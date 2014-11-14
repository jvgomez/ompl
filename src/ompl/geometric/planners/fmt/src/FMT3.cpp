/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2013, Autonomous Systems Laboratory, Stanford University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of Stanford University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Authors: Ashley Clark (Stanford) and Wolfgang Pointner (AIT) */
/* Co-developers: Brice Rebsamen (Stanford), Tim Wheeler (Stanford)
                  Edward Schmerling (Stanford), and Javier V. Gómez (UC3M - Stanford)*/
/* Algorithm design: Lucas Janson (Stanford) and Marco Pavone (Stanford) */
/* Acknowledgements for insightful comments: Oren Salzman (Tel Aviv University),
 *                                           Joseph Starek (Stanford) */

#include <limits>
#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <ompl/datastructures/BinaryHeap.h>
#include <ompl/tools/config/SelfConfig.h>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/geometric/planners/fmt/FMT3.h>

#include <fstream>
#include <ompl/base/spaces/RealVectorStateSpace.h>

ompl::geometric::FMT3::FMT3(const base::SpaceInformationPtr &si)
    : base::Planner(si, "FMT3")
    , numSamples_(1000)
    , collisionChecks_(0)
    , nearestK_(false)
    , cacheCC_(true)
    , heuristics_(false)
    , radiusMultiplier_(1.1)
    , progressiveFMT_(true)
{
    // An upper bound on the free space volume is the total space volume; the free fraction is estimated in sampleFree
    freeSpaceVolume_ = si_->getStateSpace()->getMeasure();
    lastGoalMotion_ = NULL;

    specs_.approximateSolutions = false;
    specs_.directed = false;

    ompl::base::Planner::declareParam<unsigned int>("num_samples", this, &FMT3::setNumSamples, &FMT3::getNumSamples, "10:10:1000000");
    ompl::base::Planner::declareParam<double>("radius_multiplier", this, &FMT3::setRadiusMultiplier, &FMT3::getRadiusMultiplier, "0.1:0.05:50.");
    ompl::base::Planner::declareParam<bool>("nearest_k", this, &FMT3::setNearestK, &FMT3::getNearestK, "0,1");
    ompl::base::Planner::declareParam<bool>("cache_cc", this, &FMT3::setCacheCC, &FMT3::getCacheCC, "0,1");
    ompl::base::Planner::declareParam<bool>("heuristics", this, &FMT3::setHeuristics, &FMT3::getHeuristics, "0,1");
}

ompl::geometric::FMT3::~FMT3()
{
    freeMemory();
}

void ompl::geometric::FMT3::setup()
{
    Planner::setup();

    /* Setup the optimization objective. If no optimization objective was
       specified, then default to optimizing path length as computed by the
       distance() function in the state space */
    if (pdef_->hasOptimizationObjective())
        opt_ = pdef_->getOptimizationObjective();
    else
    {
        OMPL_INFORM("%s: No optimization objective specified. Defaulting to optimizing path length.", getName().c_str());
        opt_.reset(new base::PathLengthOptimizationObjective(si_));
    }
    Open_.getComparisonOperator().opt_ = opt_.get();
    Open_.getComparisonOperator().heuristics_ = heuristics_;

    if (!nn_)
        nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion*>(si_->getStateSpace()));
    nn_->setDistanceFunction(boost::bind(&FMT3::distanceFunction, this, _1, _2));
}

void ompl::geometric::FMT3::freeMemory()
{
    if (nn_)
    {
        std::vector<Motion*> motions;
        motions.reserve(nn_->size());
        nn_->list(motions);
        for (unsigned int i = 0 ; i < motions.size() ; ++i)
        {
            si_->freeState(motions[i]->getState());
            delete motions[i];
        }
    }
}

void ompl::geometric::FMT3::clear()
{
    Planner::clear();
    lastGoalMotion_ = NULL;
    sampler_.reset();
    freeMemory();
    if (nn_)
        nn_->clear();
    Open_.clear();
    neighborhoods_.clear();

    collisionChecks_ = 0;
}

void ompl::geometric::FMT3::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);
    std::vector<Motion*> motions;
    nn_->list(motions);

    if (lastGoalMotion_)
        data.addGoalVertex(base::PlannerDataVertex(lastGoalMotion_->getState()));

    unsigned int size = motions.size();
    for (unsigned int i = 0; i < size; ++i)
    {
        if (motions[i]->getParent() == NULL)
            data.addStartVertex(base::PlannerDataVertex(motions[i]->getState()));
        else
            data.addEdge(base::PlannerDataVertex(motions[i]->getParent()->getState()),
                         base::PlannerDataVertex(motions[i]->getState()));
    }
}

void ompl::geometric::FMT3::saveNeighborhood(Motion *m)
{
    // Check to see if neighborhood has not been saved yet
    if (neighborhoods_.find(m) == neighborhoods_.end())
    {
        std::vector<Motion*> nbh;
        if (nearestK_)
            nn_->nearestK(m, NNk_, nbh);
        else
            nn_->nearestR(m, NNr_, nbh);
        if (!nbh.empty())
        {
            // Save the neighborhood but skip the first element, since it will be motion m
            neighborhoods_[m] = std::vector<Motion*>(nbh.size() - 1, 0);
            std::copy(nbh.begin() + 1, nbh.end(), neighborhoods_[m].begin());
        }
        else
        {
            // Save an empty neighborhood
            neighborhoods_[m] = std::vector<Motion*>(0);
        }
    } // If neighborhood hadn't been saved yet
}

// Calculate the unit ball volume for a given dimension
double ompl::geometric::FMT3::calculateUnitBallVolume(const unsigned int dimension) const
{
    if (dimension == 0)
        return 1.0;
    else if (dimension == 1)
        return 2.0;
    return 2.0 * boost::math::constants::pi<double>() / dimension
            * calculateUnitBallVolume(dimension - 2);
}

double ompl::geometric::FMT3::calculateRadius(const unsigned int dimension, const unsigned int n) const
{
    double a = 1.0 / (double)dimension;
    double unitBallVolume = calculateUnitBallVolume(dimension);

    return radiusMultiplier_ * 2.0 * std::pow(a, a) * std::pow(freeSpaceVolume_ / unitBallVolume, a) * std::pow(log((double)n) / (double)n, a);
}

void ompl::geometric::FMT3::sampleFree(const base::PlannerTerminationCondition &ptc)
{
    unsigned int nodeCount = 0;
    unsigned int sampleAttempts = 0;
    Motion *motion = new Motion(si_);

    // Sample numSamples_ number of nodes from the free configuration space
    while (nodeCount < numSamples_ && !ptc)
    {
        sampler_->sampleUniform(motion->getState());
        sampleAttempts++;

        bool collision_free = si_->isValid(motion->getState());

        if (collision_free)
        {
            nodeCount++;
            nn_->add(motion);
            motion = new Motion(si_);
        } // If collision free
    } // While nodeCount < numSamples
    si_->freeState(motion->getState());
    delete motion;

    // 95% confidence limit for an upper bound for the true free space volume
    freeSpaceVolume_ = boost::math::binomial_distribution<>::find_upper_bound_on_p(sampleAttempts, nodeCount, 0.05) * si_->getStateSpace()->getMeasure();
}

void ompl::geometric::FMT3::assureGoalIsSampled(const ompl::base::GoalSampleableRegion *goal)
{
    // Ensure that there is at least one node near each goal
    while (const base::State *goalState = pis_.nextGoal())
    {
        Motion *gMotion = new Motion(si_);
        si_->copyState(gMotion->getState(), goalState);

        std::vector<Motion*> nearGoal;
        nn_->nearestR(gMotion, goal->getThreshold(), nearGoal);

        // If there is no node in the goal region, insert one
        if (nearGoal.empty())
        {
            OMPL_DEBUG("No state inside goal region");
            if (si_->getStateValidityChecker()->isValid(gMotion->getState()))
            {
                nn_->add(gMotion);
                goalState_ = gMotion->getState();
            }
            else
            {
                si_->freeState(gMotion->getState());
                delete gMotion;
            }
        }
        else // There is already a sample in the goal region
        {
            goalState_ = nearGoal[0]->getState();
            si_->freeState(gMotion->getState());
            delete gMotion;
        }
    } // For each goal
}

ompl::base::PlannerStatus ompl::geometric::FMT3::solve(const base::PlannerTerminationCondition &ptc)
{
    if (lastGoalMotion_) {
        OMPL_INFORM("solve() called before clear(); returning previous solution");
        traceSolutionPathThroughTree(lastGoalMotion_);
        OMPL_DEBUG("Final path cost: %f", lastGoalMotion_->getCost().value());
        return base::PlannerStatus(true, false);
    }
    else if (Open_.size() > 0)
    {
        OMPL_INFORM("solve() called before clear(); no previous solution so starting afresh");
        clear();
    }

    checkValidity();
    base::GoalSampleableRegion *goal = dynamic_cast<base::GoalSampleableRegion*>(pdef_->getGoal().get());
    Motion *initMotion = NULL;

    if (!goal)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }

    // Add start states to V (nn_) and Open
    while (const base::State *st = pis_.nextStart())
    {
        initMotion = new Motion(si_);
        si_->copyState(initMotion->getState(), st);
        Open_.insert(initMotion);
        initMotion->setSetType(Motion::SET_OPEN);
        initMotion->setCost(opt_->initialCost(initMotion->getState()));
        nn_->add(initMotion); // V <-- {x_init}
    }

    if (!initMotion)
    {
        OMPL_ERROR("Start state undefined");
        return base::PlannerStatus::INVALID_START;
    }

    // Sample N free states in the configuration space
    if (!sampler_)
        sampler_ = si_->allocStateSampler();
    sampleFree(ptc);
    assureGoalIsSampled(goal);
    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());

    // Calculate the nearest neighbor search radius
    if (nearestK_)
    {
        NNk_ = std::ceil(std::pow(2.0 * radiusMultiplier_, (double)si_->getStateDimension()) *
                        (boost::math::constants::e<double>() / (double)si_->getStateDimension()) *
                        log((double)nn_->size()));
        OMPL_DEBUG("Using nearest-neighbors k of %d", NNk_);
    }
    else
    {
        NNr_ = calculateRadius(si_->getStateDimension(), nn_->size());
        OMPL_DEBUG("Using radius of %f", NNr_);
    }

    // Execute the planner, and return early if the planner returns a failure
    bool plannerSuccess = false;
    bool successfulExpansion = false;
    Motion *z = initMotion; // z <-- xinit
    saveNeighborhood(z);

    while (!ptc)
    {
        if ((plannerSuccess = goal->isSatisfied(z->getState())))
            break;

        successfulExpansion = expandTreeFromNode(z);

        if (!progressiveFMT_ && !successfulExpansion)
            break;
        else if (progressiveFMT_ && !successfulExpansion)
        {
            // Sample and connect samples to tree only if there is
            // a possibility to connect to unvisited nodes.
            std::vector<Motion*>       nbh;
            std::vector<base::Cost>    costs;
            std::vector<base::Cost>    incCosts;
            std::vector<std::size_t>   sortedCostIndices;

            // our functor for sorting nearest neighbors
            CostIndexCompare compareFn(costs, *opt_);

            Motion *m = new Motion(si_);
            while (!ptc && Open_.empty())
            {
                sampler_->sampleUniform(m->getState());
                // TODO: benchmark with and without isValid.
                if(!si_->isValid(m->getState()))
                    continue;

                // Does the new sample has a unvisited node as neighbor?
                nn_->nearestR(m, NNr_, nbh);
                bool connects = false;
                for(std::size_t j = 0; j < nbh.size(); ++j)
                {
                    if(nbh[j]->getSetType() == Motion::SET_UNVISITED)
                    {
                        connects = true;
                        break;
                    }
                }

                std::vector<Motion*> yNear;
                if (connects)
                {
                    // Get neighbours in the tree.
                    yNear.reserve(nbh.size());
                    for(std::size_t j = 0; j < nbh.size(); ++j)
                    {
                        if(nbh[j]->getSetType() == Motion::SET_CLOSED)
                            yNear.push_back(nbh[j]);
                    }

                    // cache for distance computations
                    //
                    // Our cost caches only increase in size, so they're only
                    // resized if they can't fit the current neighborhood
                    if (costs.size() < yNear.size())
                    {
                        costs.resize(yNear.size());
                        incCosts.resize(yNear.size());
                        sortedCostIndices.resize(yNear.size());
                    }

                    // Finding the nearest neighbor to connect to
                    // By default, neighborhood states are sorted by cost, and collision checking
                    // is performed in increasing order of cost
                    //
                    // calculate all costs and distances
                    for (std::size_t i = 0 ; i < yNear.size(); ++i)
                    {
                        incCosts[i] = opt_->motionCost(yNear[i]->getState(), m->getState());
                        costs[i] = opt_->combineCosts(yNear[i]->getCost(), incCosts[i]);
                    }

                    // sort the nodes
                    //
                    // we're using index-value pairs so that we can get at
                    // original, unsorted indices
                    for (std::size_t i = 0; i < yNear.size(); ++i)
                        sortedCostIndices[i] = i;
                    std::sort(sortedCostIndices.begin(), sortedCostIndices.begin() + yNear.size(),
                              compareFn);

                    // collision check until a valid motion is found
                   for (std::vector<std::size_t>::const_iterator i = sortedCostIndices.begin();
                        i != sortedCostIndices.begin() + yNear.size();
                        ++i)
                   {
                       if(si_->checkMotion(yNear[*i]->getState(), m->getState()))
                       {
                           m->setParent(yNear[*i]);
                           yNear[*i]->children.push_back(m);
                           const base::Cost incCost = opt_->motionCost(yNear[*i]->getState(), m->getState());
                           m->setCost(opt_->combineCosts(yNear[*i]->getCost(), incCost));
                           m->setHeuristicCost(opt_->motionCostHeuristic(m->getState(), goalState_));
                           m->setSetType(Motion::SET_OPEN);

                           nn_->add(m);
                           saveNeighborhood(m);
                           updateNeighborhood(m, nbh, NNr_);
                           Open_.insert(m);
                           z = m;
                           break;
                       }
                   }

                } // if connects
            } // while (!ptc && Open_.empty())
        } // else if (progressiveFMT_ && !successfulExpansion)
    } // While not at goal

    if (plannerSuccess)
    {
        lastGoalMotion_ = z;
        OMPL_DEBUG("Final path cost: %f", lastGoalMotion_->getCost().value());
        traceSolutionPathThroughTree(lastGoalMotion_);
        return base::PlannerStatus(true, false);
    } // if plannerSuccess
    else
    {
        // Planner terminated without accomplishing goal
        // Execute second form of anytime here.
        return base::PlannerStatus(false, false);
    }
}

void ompl::geometric::FMT3::traceSolutionPathThroughTree(Motion *goalMotion)
{
    std::vector<Motion*> mpath;
    Motion *solution = goalMotion;

    // Construct the solution path
    while (solution != NULL)
    {
        mpath.push_back(solution);
        solution = solution->getParent();
    }

    // Set the solution path
    PathGeometric *path = new PathGeometric(si_);
    int mPathSize = mpath.size();
    for (int i = mPathSize - 1 ; i >= 0 ; --i)
        path->append(mpath[i]->getState());
    pdef_->addSolutionPath(base::PathPtr(path), false, -1.0, getName());
}

bool ompl::geometric::FMT3::expandTreeFromNode(Motion *&z)
{
    // Find all nodes that are near z, and also in set Unvisited

    std::vector<Motion*> xNear;
    const std::vector<Motion*> &zNeighborhood = neighborhoods_[z];
    const unsigned int zNeighborhoodSize = zNeighborhood.size();
    xNear.reserve(zNeighborhoodSize);

    for (unsigned int i = 0; i < zNeighborhoodSize; ++i)
    {
        Motion *x = zNeighborhood[i];
        if (x->getSetType() == Motion::SET_UNVISITED)
        {
            saveNeighborhood(x);
            if (nearestK_)
            {
                // Only include neighbors that are mutually k-nearest
                // Relies on NN datastructure returning k-nearest in sorted order
                const base::Cost connCost = opt_->motionCost(z->getState(), x->getState());
                const base::Cost worstCost = opt_->motionCost(neighborhoods_[x].back()->getState(), x->getState());

                if (opt_->isCostBetterThan(worstCost, connCost))
                    continue;
                else
                    xNear.push_back(zNeighborhood[i]);
            }
            else
                xNear.push_back(x);
        }
    }

    // For each node near z and in set Unvisited, attempt to connect it to set Open
    std::vector<Motion*> yNear;
    const unsigned int xNearSize = xNear.size();
    for (unsigned int i = 0 ; i < xNearSize; ++i)
    {
        Motion *x = xNear[i];

        // Find all nodes that are near x and in set Open
        const std::vector<Motion*> &xNeighborhood = neighborhoods_[x];

        const unsigned int xNeighborhoodSize = xNeighborhood.size();
        yNear.reserve(xNeighborhoodSize);
        for (unsigned int j = 0; j < xNeighborhoodSize; ++j)
        {
            if (xNeighborhood[j]->getSetType() == Motion::SET_OPEN)
                yNear.push_back(xNeighborhood[j]);
        }

        // Find the lowest cost-to-come connection from Open to x
        base::Cost cMin(std::numeric_limits<double>::infinity());
        Motion *yMin = getBestParent(x, yNear, cMin);
        yNear.clear();

        // If an optimal connection from Open to x was found
        if (yMin != NULL)
        {
            bool collision_free = false;
            if (cacheCC_)
            {
                if (!yMin->alreadyCC(x))
                {
                    collision_free = si_->checkMotion(yMin->getState(), x->getState());
                    ++collisionChecks_;
                    // Due to FMT3* design, it is only necessary to save unsuccesful
                    // connection attemps because of collision
                    if (!collision_free)
                        yMin->addCC(x);
                }
            }
            else
            {
                ++collisionChecks_;
                collision_free = si_->checkMotion(yMin->getState(), x->getState());
            }

            if (collision_free)
            {
                // Add edge from yMin to x
                x->setParent(yMin);
                x->setCost(cMin);
                x->setHeuristicCost(opt_->motionCostHeuristic(x->getState(), goalState_));
                yMin->children.push_back(x);
                // Add x to Open
                Open_.insert(x);
                // Remove x from Unvisited and add to Open
                x->setSetType(Motion::SET_OPEN);
            }
        } // An optimal connection from Open to x was found
    } // For each node near z and in set Unvisited, try to connect it to set Open

    // Update Open
    Open_.pop();
    z->setSetType(Motion::SET_CLOSED);

    if (Open_.empty())
    {
        OMPL_INFORM("Open is empty before path was found --> no feasible path exists");
        return false;
    }

    // Take the top of Open as the new z
    z = Open_.top()->data;

    return true;
}

std::string ompl::geometric::FMT3::getCollisionCheckCount() const
{
    return boost::lexical_cast<std::string>(collisionChecks_);
}

ompl::geometric::FMT3::Motion* ompl::geometric::FMT3::getBestParent(Motion *m, std::vector<Motion*> &neighbors, base::Cost &cMin)
{
    Motion *min = NULL;
    const unsigned int neighborsSize = neighbors.size();
    for (unsigned int j = 0; j < neighborsSize; ++j)
    {
        const base::State *s = neighbors[j]->getState();
        const base::Cost dist = opt_->motionCost(s, m->getState());
        const base::Cost cNew = opt_->combineCosts(neighbors[j]->getCost(), dist);

        if (opt_->isCostBetterThan(cNew, cMin))
        {
            min = neighbors[j];
            cMin = cNew;
        }
    }
    return min;
}

void ompl::geometric::FMT3::findLeafNodes(std::vector<Motion*> &leaves)
{
    // TODO: can be optimized by adding a flag to the motion.
    // TODO: it may also be faster to check by childrens (tree traversal) instead of listing nodes.
    std::vector<Motion*> tree;
    nn_->list(tree);

    for(std::size_t i = 0; i < tree.size(); ++i)
    {
        if(tree[i]->getSetType() == Motion::SET_CLOSED && tree[i]->children.size() == 0)
            leaves.push_back(tree[i]);
    }
}

// TODO: this could be better implemented by returning something in saveNeighborhood()
// TODO: note that the new neighbours added are not sorted.
void ompl::geometric::FMT3::updateNeighborhood(Motion *m, const std::vector<Motion*> nbh, const double r)
{
    for(std::size_t i = 0; i < nbh.size(); ++i)
    {
        // If the neihborhood has not been saved yet, m will be automatically added
        // since it is already included in the NN datastructure.
        if (neighborhoods_.find(nbh[i]) == neighborhoods_.end())
        {
            std::vector<Motion*> nbh2;
            nn_->nearestR(nbh[i], r, nbh2);
            if (!nbh2.empty())
            {
                // Save the neighborhood but skip the first element, since it will be motion m
                neighborhoods_[nbh[i]] = std::vector<Motion*>(nbh2.size() - 1, 0);
                std::copy(nbh2.begin() + 1, nbh2.end(), neighborhoods_[nbh[i]].begin());
            }
            else
            {
                // Save an empty neighborhood
                neighborhoods_[nbh[i]] = std::vector<Motion*>(0);
            }
        }
        // Neighborhood already saved, so we have to add m to it.
        // TODO : it should be included sorted!
        else
        {
            neighborhoods_[nbh[i]].push_back(m);
        }
    }
}


void ompl::geometric::FMT3::saveTree(const std::string &filename)
{
    std::ofstream ofs;
    ofs.open(filename.c_str(), std::ofstream::trunc);
    ofs << std::setprecision(6);
    std::vector<Motion *> motions;
    nn_->list(motions);

    for (size_t i = 0; i < motions.size(); ++i)
    {
        if (motions[i]->getParent())
            ofs << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[0] << "\t"
                << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[1] << "\t"
                << motions[i]->getParent()->getState()->as<base::RealVectorStateSpace::StateType>()->values[0] << "\t"
                << motions[i]->getParent()->getState()->as<base::RealVectorStateSpace::StateType>()->values[1] << "\t"
                << motions[i]->getCost() << "\t" <<

                   std::endl;
                   /*motions[i]->getHeuristicCost() << "\t"
                << opt_->combineCosts(motions[i]->getCost(), motions[i]->getHeuristicCost()) << std::endl;*/
        else
            ofs << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[0] << "\t"
                << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[1] << "\t"
                << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[0] << "\t"
                << motions[i]->getState()->as<base::RealVectorStateSpace::StateType>()->values[1] << "\t"
                << 0 <<
                std::endl;
                //"\t" << 0 << "\t"<< 0 << std::endl;
    }

    ofs.close();
}
