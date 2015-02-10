
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <ompl/datastructures/BinaryHeap.h>
#include <ompl/tools/config/SelfConfig.h>

#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <ompl/base/objectives/PathLengthOptimizationObjective.h>
#include <ompl/geometric/planners/fmt/BFMT.h>
// #include "BiDirectionalFMT.h"

namespace ompl {
namespace geometric {

BFMT::BFMT(const base::SpaceInformationPtr &si)
    : base::Planner(si, "BFMT")
    , numSamples_(1000)
    , radiusMultiplier_(1.0)
    , freeSpaceVolume_(si_->getStateSpace()->getMeasure())      // An upper bound on the free space volume is the total space volume; the free fraction is estimated in sampleFree
    , collisionChecks_(0)
    , nearestK_(true)
    , NNr(0)
    , NNk(0)
    , tree_(FWD)
    , exploration_(SWAP_EVERY_TIME)
    , termination_(OPTIMALITY)
    , precomputeNN_(false)
{
    specs_.approximateSolutions = false;
    specs_.directed             = false;
    
    ompl::base::Planner::declareParam<unsigned int>("num_samples", this, &BFMT::setNumSamples, &BFMT::getNumSamples, "10:10:1000000" );
    ompl::base::Planner::declareParam<double>("radius_multiplier", this, &BFMT::setRadiusMultiplier, &BFMT::getRadiusMultiplier, "0.1:0.05:50." );
    ompl::base::Planner::declareParam<bool>("nearest_k", this, &BFMT::setNearestK, &BFMT::getNearestK, "0,1" );
    ompl::base::Planner::declareParam<bool>("balanced", this, &BFMT::setExploration, &BFMT::getExploration, "0,1" );
    ompl::base::Planner::declareParam<bool>("optimality", this, &BFMT::setTermination, &BFMT::getTermination, "0,1" );
}
                
BFMT::BFMT(const base::SpaceInformationPtr &si, const bool precompute_NN)
    : base::Planner(si, "BFMT")
    , numSamples_(1000)
    , radiusMultiplier_(1.0)
    , freeSpaceVolume_(si_->getStateSpace()->getMeasure())      // An upper bound on the free space volume is the total space volume; the free fraction is estimated in sampleFree
    , collisionChecks_(0)
    , nearestK_(true)
    , NNr(0)
    , NNk(0)
    , tree_(FWD)
    , exploration_(SWAP_EVERY_TIME)
    , termination_(OPTIMALITY)
    , precomputeNN_(precompute_NN)
{
}

BFMT::BFMT(const base::SpaceInformationPtr &si, const enum ExploreType exploration, 
        const enum TerminateType termination, const bool precompute_NN)
    : base::Planner(si, "BFMT")
    , numSamples_(1000)
    , radiusMultiplier_(1.0)
    , freeSpaceVolume_(si_->getStateSpace()->getMeasure())      // An upper bound on the free space volume is the total space volume; the free fraction is estimated in sampleFree
    , collisionChecks_(0)
    , nearestK_(true)
    , NNr(0)
    , NNk(0)
    , tree_(FWD)
    , exploration_(exploration)
    , termination_(termination)
    , precomputeNN_(precompute_NN)
{
}

ompl::geometric::BFMT::~BFMT() {
    freeMemory();
}

void BFMT::setup(void) {
    ompl::base::Planner::setup();

    /* Setup the optimization objective. If no optimization objective was
       specified, then default to optimizing path length as computed by the
       distance() function in the state space */
    if (pdef_->hasOptimizationObjective()) {
        opt_ = pdef_->getOptimizationObjective();
    } else {
        OMPL_INFORM("%s: No optimization objective specified. Defaulting to optimizing path length.", getName().c_str());
        opt_.reset(new base::PathLengthOptimizationObjective(si_));
    }
    // H.getComparisonOperator().opt_ = opt_.get();
    
    if (!nn_) {
        nn_.reset(new NearestNeighborsGNAT<BiDirMotion*>());
    }
    nn_->setDistanceFunction(boost::bind(&BFMT::distanceFunction, this, _1, _2));
}
    
void BFMT::freeMemory() {
    if (nn_) {
        BiDirMotionPtrs motions;
        nn_->list(motions);
        for (unsigned int i = 0 ; i < motions.size() ; ++i) {
            si_->freeState(motions[i]->getState());
            delete motions[i];
        }
    }
}

void BFMT::clear() {
    Planner::clear();
    sampler_.reset();
    freeMemory();
    if (nn_)
        nn_->clear();
    H[FWD].clear();
    H[REV].clear();
    H_elements[FWD].clear();
    H_elements[REV].clear();
    neighborhoods_.clear();
    collisionChecks_ = 0;
}

void BFMT::getPlannerData( base::PlannerData &data ) const {
    base::Planner::getPlannerData(data);
    BiDirMotionPtrs motions;
    nn_->list(motions);

    int numStartNodes       = 0;
    int numGoalNodes        = 0;
    int numEdges            = 0;
    int numFwdEdges         = 0;
    int numRevEdges         = 0;

    int fwd_tree_tag        = 1;
    int rev_tree_tag        = 2;

    for ( unsigned int k = 0; k < motions.size(); ++k ) {
        BiDirMotion* motion = motions[k];
        bool inFwdTree      = ( motion->currentSet_[FWD] != BiDirMotion::SET_W ); //|| Nfwd_children > 0 );

        // For samples added to the fwd tree, add incoming edges (from fwd tree parent)
        if ( inFwdTree ) {
            if ( motion->parent_[FWD] == NULL ) {
                // Motion is a forward tree root node
                ++numStartNodes;
            } else {
                bool success = data.addEdge(
                    base::PlannerDataVertex(motion->parent_[FWD]->getState(), fwd_tree_tag),
                    base::PlannerDataVertex(motion->getState(), fwd_tree_tag)
                );
                if (success) {
                    ++numFwdEdges;
                    ++numEdges;
                }
            }
        }
    }

    // The edges in the goal tree are reversed so that they are in the same direction as start tree
    for ( unsigned int k = 0; k < motions.size(); ++k ) {
        BiDirMotion* motion = motions[k];
        bool inRevTree      = ( motion->currentSet_[REV] != BiDirMotion::SET_W ); //|| Nrev_children > 0 );

        // For samples added to a tree, add incoming edges (from fwd tree parent)
        if ( inRevTree ) {
            if ( motion->parent_[REV] == NULL ) {
                // Motion is a reverse tree root node
                ++numGoalNodes;
            } else {
                bool success = data.addEdge(
                    base::PlannerDataVertex(motion->getState(), rev_tree_tag),
                    base::PlannerDataVertex(motion->parent_[REV]->getState(), rev_tree_tag)
                );
                if (success) {
                    ++numRevEdges;
                    ++numEdges;
                }
            }
        }
    }
}

void BFMT::saveNeighborhood( boost::shared_ptr< NearestNeighbors<BiDirMotion*> > nn, 
        BiDirMotion* m) {
    
    // Check if neighborhood has already been saved
    if ( neighborhoods_.find(m) == neighborhoods_.end() ) {
        BiDirMotionPtrs neighborhood;
        if (nearestK_) {
            nn_->nearestK(m, NNk, neighborhood);
        } else {
            nn_->nearestR(m, NNr, neighborhood);
        }
        
        if (!neighborhood.empty()) {
            // Save the neighborhood but skip the first element (m)
            neighborhoods_[m] = std::vector<BiDirMotion*>(neighborhood.size()-1, 0);
            std::copy(neighborhood.begin()+1, neighborhood.end(), neighborhoods_[m].begin());
            
        } else {
            // Save an empty neighborhood
            neighborhoods_[m] = std::vector<BiDirMotion*>(0);
        }
    }
}


void BFMT::sampleFree( boost::shared_ptr<NearestNeighbors<BiDirMotion*> > nn,
        const base::PlannerTerminationCondition &ptc) {
    
    unsigned int nodeCount  = 0;
    BiDirMotion *motion     = new BiDirMotion( si_, &tree_ );

    while (nodeCount < numSamples_ && !ptc) {
        sampler_->sampleUniform(motion->getState());
        if (si_->isValid(motion->getState())) {     // collision checking
            ++nodeCount;
            nn->add(motion);
            motion = new BiDirMotion( si_, &tree_ );
        }
    }
    delete motion;
}

// Helper functions
double BFMT::calculateUnitBallVolume(const unsigned int dimension) const {
    if ( dimension == 0 )
        return 1.0;
    else if( dimension == 1 )
        return 2.0;
    return 2.0 * boost::math::constants::pi<double>() / dimension
            * calculateUnitBallVolume(dimension-2);
}

double BFMT::calculateRadius(const unsigned int dimension, const unsigned int n) const {
    double a = 1.0/(double)dimension;
    double lebesqueMeasure = freeSpaceVolume_;
    double unitBallVolume = calculateUnitBallVolume(dimension);

    double radius = radiusMultiplier_ * 2.0 * pow(1.0+a,a) * pow(lebesqueMeasure/unitBallVolume, a) * pow(log((double)n)/(double)n,a);
    return radius;
}

void BFMT::initializeProblem( base::GoalSampleableRegion*& goal_s) {
    checkValidity();
    if (!sampler_) {
        sampler_ = si_->allocStateSampler();
    }
    goal_s = dynamic_cast<base::GoalSampleableRegion*>( pdef_->getGoal().get() );
}


// Solve command
base::PlannerStatus BFMT::solve(const base::PlannerTerminationCondition& ptc) {

    base::GoalSampleableRegion  *goal_s;
    initializeProblem(goal_s);
    if (!goal_s) {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // All actions from this point augment the forward tree
    ///////////////////////////////////////////////////////////////////////////
    useFwdTree();
    
    // Add start states to Vfwd and Hfwd
    bool valid_initMotion = false;
    BiDirMotion *initMotion;
    while (const base::State *st = pis_.nextStart()) {
        initMotion = new BiDirMotion( si_, &tree_ );
        si_->copyState(initMotion->getState(), st);
        
        initMotion->currentSet_[REV]        = BiDirMotion::SET_W;
        nn_->add(initMotion);                                                   // S <-- {x_init}
        if (si_->isValid(initMotion->getState())) {
            // Take the first valid initial state as the forward tree root
            H_elements[FWD][initMotion]     = H[FWD].insert(initMotion);
            initMotion->currentSet_[FWD]    = BiDirMotion::SET_H;
            initMotion->cost_[FWD]          = opt_->initialCost(initMotion->getState());
            valid_initMotion                = true;
        }
    }
    
    if (!initMotion || !valid_initMotion) {
        OMPL_ERROR("Start state undefined or invalid.");
        return base::PlannerStatus::INVALID_START;
    }
    
    // Sample N free states in configuration state_
    // TODO move this up and label all N new points in Wfwd and Wrev???
    sampleFree(nn_, ptc);                                                       // S <-- SAMPLEFREE(N)
    OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());
    
    // Calculate the nearest neighbor search radius
    if (nearestK_) {
        NNk = std::ceil(std::pow(2.0 * radiusMultiplier_, (double)si_->getStateDimension()) *
                        (boost::math::constants::e<double>() / (double)si_->getStateDimension()) *
                        log((double)nn_->size()));
        OMPL_DEBUG("Using nearest-neighbors k of %d", NNk);
    } else {
        NNr = calculateRadius(si_->getStateDimension(), nn_->size());
        OMPL_DEBUG("Using radius of %f", NNr);
    }
    
    // Add goal states to Vrev and Hrev
    bool valid_goalMotion = false;
    BiDirMotion *goalMotion;
    while (const base::State *st = pis_.nextGoal()) {
        goalMotion = new BiDirMotion(si_, &tree_);
        si_->copyState(goalMotion->getState(), st);
        
        goalMotion->currentSet_[FWD]        = BiDirMotion::SET_W;
        nn_->add(goalMotion);                                                   // S <-- {x_goal}
        if (si_->isValid(goalMotion->getState())) {
            // Take the first valid goal state as the reverse tree root
            H_elements[REV][goalMotion]     = H[REV].insert(goalMotion);
            goalMotion->currentSet_[REV]    = BiDirMotion::SET_H;
            goalMotion->cost_[REV]          = opt_->terminalCost(goalMotion->getState());
            valid_goalMotion                = true;
        }
    }
    
    if (!goalMotion || !valid_goalMotion) {
        OMPL_ERROR("Goal state undefined or invalid.");
        return base::PlannerStatus::INVALID_GOAL;
    }
    
//    // Sample from the goal region and check that nearest-neighbor connections exist
//    BiDirMotion *gMotion    = new BiDirMotion( si_, &tree_ );
//    base::State* gState     = si_->allocState();
//    goal_s->sampleGoal( gMotion->getState() );
//    si_->copyState( gState, gMotion->getState() );
//    
//    BiDirMotionPtrs nearGoal;
//    nn_->nearestR(gMotion, goal_s->getThreshold(), nearGoal); // nearest neighbors
//
//    if (nearGoal.empty()) {
//        OMPL_DEBUG("No state_ inside goal region");
//        gMotion->setCurrentSet(BiDirMotion::SET_W);           // insert goal motion node into Wfwd
//    }
    
    
    ///////////////////////////////////////////////////////////////////////////
    // All actions from this point augment the reverse tree
    ///////////////////////////////////////////////////////////////////////////
    useRevTree();
    
//    // add goal states to Vrev and Hrev
//    H_elements[tree_][gMotion] = H[tree_].insert(gMotion);
//    gMotion->setCurrentSet(BiDirMotion::SET_H);
//    nn_->add(gMotion);				// V <-- {x_init}
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Plan a path
    ///////////////////////////////////////////////////////////////////////////
    BiDirMotion *connection_Point = NULL;
    bool earlyFailure = true;
    
    if (initMotion != NULL && goalMotion != NULL) {
        earlyFailure = plan( initMotion, goalMotion, connection_Point, ptc);
    } else {
        OMPL_ERROR("Initial/goal state(s) are undefined!");
    }
    
    if (earlyFailure) {
        return base::PlannerStatus(false,false);
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Save the best path (through z)
    ///////////////////////////////////////////////////////////////////////////
    if (!ptc) {
        base::Cost fwd_cost, rev_cost, connection_cost;
        
        // Construct the solution path
        useFwdTree();
        BiDirMotionPtrs path_fwd;
        tracePath(connection_Point, path_fwd);
        fwd_cost = connection_Point->getCost();

        useRevTree();
        BiDirMotionPtrs path_rev;
        tracePath(connection_Point, path_rev);
        rev_cost = connection_Point->getCost();

        // ASSUMES FROM THIS POINT THAT z = path_fwd[0] = path_rev[0]
        // Remove the first element, z, in the traced reverse path 
        // (the same as the first element in the traced forward path)
        assert( path_fwd[0] == path_rev[0] );
        assert( path_fwd[0] == connection_Point );
        if (path_rev.size() > 1) {
            connection_cost = base::Cost(rev_cost.value() - path_rev[1]->getCost().value());
            path_rev.erase(path_rev.begin());
        } else if (path_fwd.size() > 1) {
            connection_cost = base::Cost(fwd_cost.value() - path_fwd[1]->getCost().value());
            path_fwd.erase(path_fwd.begin());
        } else {
            OMPL_ERROR("Solution path traced incorrectly or otherwise constructed improperly \
                through forward/reverse trees (both paths are one node in length, each).");
        }
        
        // Adjust costs/parents in reverse tree nodes as cost/direction from forward tree root
        useFwdTree();
        path_rev[0]->setCost( base::Cost(path_fwd[0]->getCost().value() + connection_cost.value()) );       // Valid if either path_fwd[0] or path_rev[0] is erased
        path_rev[0]->setParent(path_fwd[0]);
        for ( unsigned int i = 1; i < path_rev.size(); ++i ) {
            path_rev[i]->setCost( base::Cost(fwd_cost.value() + (rev_cost.value() - path_rev[i]->getCost().value())) );   // Valid only if path_rev[0] is removed, but for loop condition prevents this from being executed
            path_rev[i]->setParent( path_rev[i-1] );
        }
        
        BiDirMotionPtrs mpath;
        std::reverse(path_rev.begin(), path_rev.end());
        mpath.reserve( path_fwd.size() + path_rev.size() );                // preallocate memory
        mpath.insert( mpath.end(), path_rev.begin(), path_rev.end() );
        mpath.insert( mpath.end(), path_fwd.begin(), path_fwd.end() );
        
        // Set the solution path
        PathGeometric *path = new PathGeometric(si_);
        for ( int i = mpath.size() - 1 ; i >= 0 ; --i ) {
            path->append(mpath[i]->getState());
        }

        static const bool    approximate                 = false;
        static const double  cost_difference_from_goal   = 0.0;
        pdef_->addSolutionPath( base::PathPtr(path), approximate, cost_difference_from_goal, getName() );
        //OMPL_DEBUG("Fwd tree path length = %i\n", path_fwd.size());
        //OMPL_DEBUG("Rev tree path length = %i\n", path_rev.size());
        //OMPL_DEBUG("z Cost-to-Come \t= ", fwd_cost.v );
        //OMPL_DEBUG("z Cost-to-Go \t= ", rev_cost.v );

        //std::cout << "fwd_cost.v: " << fwd_cost.v << std::endl;
        //std::cout << "rev_cost.v: " << rev_cost.v << std::endl;
        //std::cout << "Total path cost: " << fwd_cost.v + rev_cost.v << std::endl;
        OMPL_DEBUG("Total path cost: %f\n", fwd_cost.value() + rev_cost.value() );
        return base::PlannerStatus(true, false);
        
    } else {
        // Planner terminated without accomplishing goal
        return base::PlannerStatus(false, false);
    }
}


void BFMT::expandTreeFromNode( BiDirMotion *&z, BiDirMotion *&connection_Point ) {

    // Define Hnew and set it to NULL
    BiDirMotionPtrs        H_new;
    
    // Define Znear as all unexplored nodes in the neighborhood around z
    BiDirMotionPtrs        zNear;
    const BiDirMotionPtrs  &zNeighborhood = neighborhoods_[z];

    for ( unsigned int i = 0; i < zNeighborhood.size(); ++i ) {
        if (zNeighborhood[i]->getCurrentSet() == BiDirMotion::SET_W) {
            zNear.push_back(zNeighborhood[i]);
        }
    }
    // OMPL_DEBUG("zNbhd size = %i\n", zNeighborhood.size());
    // OMPL_DEBUG("Znear size = %i\n", zNear.size());
    
    // For each node x in Znear
    for ( unsigned int i = 0; i < zNear.size(); ++i ) {
        BiDirMotion *x = zNear.at(i);
        if (precomputeNN_ == false) {
            saveNeighborhood(nn_, x); // nearest neighbors
        }
        
        // Define Xnear as all frontier nodes in the neighborhood around the unexplored node x
        BiDirMotionPtrs xNear;
        const BiDirMotionPtrs &xNeighborhood = neighborhoods_[x];
        for ( unsigned int j = 0; j < xNeighborhood.size(); ++j ) {
            if (xNeighborhood[j]->getCurrentSet() == BiDirMotion::SET_H) {
                xNear.push_back(xNeighborhood[j]);
            }
        }

        // Find the node in Xnear with minimum cost-to-come in the current tree
        BiDirMotion* xMin   = NULL;
        double cMin         = std::numeric_limits<double>::infinity();
        for ( unsigned int j = 0; j < xNear.size(); ++j ) {

            // check if node costs are smaller than minimum
            double cNew = xNear.at(j)->getCost().value() + distanceFunction(xNear.at(j), x);

            if (cNew < cMin) {
                xMin = xNear.at(j);
                cMin = cNew;
            }
        }


        // xMin was found
        if (xMin != NULL) {
            ++collisionChecks_;
            bool collision_free = si_->checkMotion(xMin->getState(), x->getState()); // collision checking

            if (collision_free) { // motion between yMin and x is obstacle free
                // add edge from xMin to x
                x->setParent(xMin);
                x->setCost( cMin );
                xMin->getChildren().push_back(x);

                // check if new node x is in the other tree; if so, save result
                if (x->getOtherSet() != BiDirMotion::SET_W) {
                    if (connection_Point == NULL) {
                        connection_Point = x;
                        if (termination_ == FEASIBILITY) {
                            break;
                        }
                    } else {
                        if ( (connection_Point->cost_[FWD].value() + connection_Point->cost_[REV].value())
                                > (x->cost_[FWD].value() + x->cost_[REV].value()) ) {
                            connection_Point = x;
                        }
                    }
                }
                
                H_new.push_back(x);                     // add x to H_new
                //TODO remove the following line??
                x->setCurrentSet(BiDirMotion::SET_NULL);     // remove x from W
            }
        }
    } // End "for x in Znear"
    
    // Remove motion z from binary heap and map
    BiDirMotionBinHeap::Element* zElement = H_elements[tree_][z];
    H[tree_].remove(zElement);
    H_elements[tree_].erase(z);
    z->setCurrentSet(BiDirMotion::SET_NULL);
    //OMPL_DEBUG("Hnew size = %i\n", H_new.size());
    
    // add nodes in H_new to H
    for (unsigned int i = 0; i < H_new.size(); i++) {
        H_elements[tree_][H_new.at(i)] = H[tree_].insert(H_new.at(i));
        H_new.at(i)->setCurrentSet(BiDirMotion::SET_H);
    }
}


bool BFMT::plan( BiDirMotion *x_init, BiDirMotion *x_goal, 
        BiDirMotion *&connection_Point, const base::PlannerTerminationCondition& ptc) {

    // If pre-computation, find neighborhoods for all N sample nodes plus initial
    // and goal state(s).  Otherwise compute the neighborhoods of the initial and
    // goal states separately and compute the others as needed.
    BiDirMotionPtrs sampleNodes;
    nn_->list(sampleNodes);
    if (precomputeNN_ == true) {
        for (unsigned int i = 0; i < sampleNodes.size(); i++) {
            saveNeighborhood( nn_, sampleNodes[i] ); // nearest neighbors
        }
    } else {
        saveNeighborhood( nn_, x_init ); // nearest neighbors
        saveNeighborhood( nn_, x_goal ); // nearest neighbors
    }

    // Copy nodes in the sample set to Wfwd.  Overwrite the label of the initial
    // node with set H for the forward tree, since it starts in set Hfwd.
    useFwdTree();
    for (unsigned int i = 0; i < sampleNodes.size(); i++) {
        sampleNodes[i]->setCurrentSet(BiDirMotion::SET_W);
    }
    x_init->setCurrentSet(BiDirMotion::SET_H);
    
    // Copy nodes in the sample set to Wrev.  Overwrite the label of the goal
    // node with set H for the reverse tree, since it starts in set Hrev.
    useRevTree();
    for (unsigned int i = 0; i < sampleNodes.size(); i++) {
        sampleNodes[i]->setCurrentSet(BiDirMotion::SET_W);
    }
    x_goal->setCurrentSet(BiDirMotion::SET_H);

    
    ///////////////////////////////////////////////////////////////////////////
    // Expand the trees until reaching the termination condition
    ///////////////////////////////////////////////////////////////////////////
    bool earlyFailure   = false;
    bool success        = false;
    
    useFwdTree();
    BiDirMotion *z = x_init;
    
    while (success == false) {

        expandTreeFromNode( z, connection_Point );
        
        // Check if the algorithm should terminate.  Possibly redefines connection_Point.
        if ( termination(z, connection_Point, ptc) ) {
            success = true;
        } else {
            if (H[tree_].empty() && H[(tree_+1) % 2].empty()) {
                OMPL_INFORM("Both H are empty before path was found --> no feasible path exists");
                earlyFailure = true;
                return earlyFailure;
            }

            chooseTreeAndExpansionNode(z);
        }
    }

    earlyFailure = false;
    return earlyFailure;
}

bool BFMT::termination( BiDirMotion *&z, BiDirMotion *&connection_Point, const base::PlannerTerminationCondition& ptc ) {
    bool terminate = false;
    switch (termination_) {
        case FEASIBILITY:
            // Test if a connection point was found during tree expansion
            return (connection_Point != NULL || ptc);
            break;

        case OPTIMALITY:
            // Test if z is in SET_NULL (interior) of other tree
            if (ptc) {
                terminate           = true;
            } else if (z->getOtherSet() == BiDirMotion::SET_NULL) {
                terminate           = true;
            }
            break;
    };
    return terminate;
}

// Choose exploration tree and node z to expand
void BFMT::chooseTreeAndExpansionNode( BiDirMotion *&z ) {
    BiDirMotion *z1, *z2;

    //std::vector<std::string> treestrings;
    //treestrings.push_back("FWD");
    //treestrings.push_back("REV");

    switch (exploration_) {
        case SWAP_EVERY_TIME:
            //std::cout << "Swap! ";
            if (H[(tree_+1) % 2].empty()) {
                //std::cout << "(Tree " << treestrings[(tree_+1) % 2] << " empty!!)";
                z   = H[tree_].top()->data;         //*< Continue expanding the current tree (not empty by exit condition in plan())
            } else {
                //std::cout << "Old Tree = " << treestrings[tree_];
                z   = H[(tree_+1) % 2].top()->data; //*< Take top of opposite tree heap as new z
                swapTrees();                        //*< Swap to the opposite tree
                //std::cout << ", New Tree = " << treestrings[tree_];
            }
            break;

        case CHOOSE_SMALLEST_Z:
            //std::cout << "Balance! ";
            if (H[(tree_+1) % 2].empty()) {
                //std::cout << "(Tree " << treestrings[(tree_+1) % 2] << " empty!!)";
                z   = H[tree_].top()->data;         //*< Continue expanding the current tree (not empty by exit condition in plan())
            } else if (H[tree_].empty()) {
                //std::cout << "(Tree " << treestrings[tree_] << " empty!!)";
                z = H[(tree_+1) % 2].top()->data;   //*< Take top of opposite tree heap as new z
                swapTrees();                        //*< Swap to the opposite tree
            } else {
                z1  = H[tree_].top()->data;
                z2  = H[(tree_+1) % 2].top()->data;

                //std::cout << treestrings[tree_] << " Tree Cost = " << z1->getCost() << ", ";
                //std::cout << treestrings[(tree_+1) % 2] << " Tree Cost = " << z2->getOtherCost();
                //std::cout << "Choosing tree to expand..." << std::endl;
                if ( z1->getCost().value() < z2->getOtherCost().value() ) {
                    z = z1;
                } else {
                    z = z2;
                    swapTrees();
                }
            }
            break;
    };
    //std::cout << ", Current Tree = " << treestrings[tree_] << ", Cost = " << z->getCost() << std::endl;
}

// Trace a path of nodes along a tree towards the root (forward or reverse)
void BFMT::tracePath( BiDirMotion *z, BiDirMotionPtrs& path ) {
    BiDirMotion* solution = z;

    while (solution != NULL) {
        path.push_back(solution);
        solution = solution->getParent();
    }
}

void BFMT::swapTrees() {
    tree_ = (TreeType) ( (((int) tree_) + 1) % 2 );
}

std::string BFMT::getCollisionCheckCount() const {
    return boost::lexical_cast<std::string>(collisionChecks_);
}

std::string BFMT::getNodeCount() const {
    return getExploredNodeCount();
}

std::string BFMT::getExploredNodeCount() const {
    // THIS IS VALID FOR "GETNODECOUNT" BUT FORGETS TO INCLUDE NODES IN SET_W 
    // THAT WERE EXAMINED IN "EXPANDTREEFROMNODE" DURING THE SEARCH OVER
    // NEAREST NEIGHBORS; SOME NODES IN SET_W IN THE NN OF FRONTIER NODES ARE MISSING
    unsigned int Nexplored = 0;
    BiDirMotionPtrs nodes;
    nn_->list(nodes);
    for ( unsigned int i = 0; i < nodes.size(); ++i ) {
	if (nodes[i]->currentSet_[FWD] != BiDirMotion::SET_W || nodes[i]->currentSet_[REV] != BiDirMotion::SET_W) {
	    ++Nexplored;
        }
    }
    return boost::lexical_cast<std::string>(Nexplored);
}

}   // End "geometric" namespace
}   // End "ompl" namespace
