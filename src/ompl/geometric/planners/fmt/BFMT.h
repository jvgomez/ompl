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

/* Authors: Joseph Starek (Stanford) */
/* Co-developers: Javier V Gomez (UC3M)*/
/* Algorithm design: Joseph Starek (Stanford), Ed Schmerling (Stanford), Lucas Janson (Stanford) and Marco Pavone (Stanford) */
/* Acknowledgements for insightful comments: Ashley Clark (Stanford) */

#ifndef BIDIRECTIONALFMT_H
#define BIDIRECTIONALFMT_H

#include <ompl/geometric/planners/PlannerIncludes.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include <ompl/datastructures/BinaryHeap.h>
#include <ompl/base/OptimizationObjective.h>
#include <map>

namespace ompl {
    namespace geometric {
        
        class BFMT : public ompl::base::Planner {
        public:
            enum TreeType { FWD = 0, REV = 1 };                                 //*< Tree identifier
            enum ExploreType { SWAP_EVERY_TIME = 0, CHOOSE_SMALLEST_Z = 1 };    //*< Exploration strategy identifier
            enum TerminateType { FEASIBILITY = 0, OPTIMALITY = 2 };             //*< Termination strategy identifier
            
            BFMT(const base::SpaceInformationPtr &si);
                
            BFMT(const base::SpaceInformationPtr &si, const bool precompute_NN);
            
            BFMT(const base::SpaceInformationPtr &si, const enum ExploreType exploration, 
                    const enum TerminateType termination, const bool precompute_NN);
            
            virtual ~BFMT();
            
            virtual void setup(void);
            
            virtual base::PlannerStatus solve(const base::PlannerTerminationCondition& ptc);
            
            virtual void clear();
            
            virtual void getPlannerData( base::PlannerData &data ) const;
            
            // Data collection functions (for modified OMPL.app GUI)
            void setNumSamples(const unsigned int numSamples) {
                numSamples_ = numSamples;
            }

            unsigned int getNumSamples() const {
                return numSamples_;
            }

            std::string getCollisionCheckCount() const;
            std::string getNodeCount() const;
            std::string getExploredNodeCount() const;
            // End data collection functions

            // Other helper functions
            void setRadiusMultiplier(const double radiusMultiplier) {
                if (radiusMultiplier <= 0.0)
                    throw Exception("Radius multiplier must be greater than zero");
                radiusMultiplier_ = radiusMultiplier;
            }

            double getRadiusMultiplier() const {
                return radiusMultiplier_;
            }

            void setNearestK(bool nearestK) {
                nearestK_ = nearestK;
            }

            bool getNearestK() const {
                return nearestK_;
            }

            void setExploration(bool balanced){
                exploration_ = SWAP_EVERY_TIME;
                if (balanced) {
                    exploration_ = CHOOSE_SMALLEST_Z;
                }
            }
            
            bool getExploration(void) const {
                return (exploration_ == CHOOSE_SMALLEST_Z);
            }
            
            void setTermination(bool optimality) {
                termination_ = FEASIBILITY;
                if (optimality) {
                    termination_ = OPTIMALITY;
                }
            }
            
            bool getTermination(void) const {
                return (termination_ == OPTIMALITY);
            }
            
            void setFreeSpaceVolume(const double freeSpaceVolume) {
                if (freeSpaceVolume < 0.0)
                    throw Exception("Free space volume should be greater than zero");
                freeSpaceVolume_ = freeSpaceVolume;
            }

            double getFreeSpaceVolume() const {
                return freeSpaceVolume_;
            }

            /** \brief Activates the cost to go heuristics when ordering the heap */
            void setHeuristics (bool h)
            {
                heuristics_ = h;
            }

            /** \brief Returns true if the heap is ordered taking into account
            cost to go heuristics */
            bool getHeuristics() const
            {
                return heuristics_;
            }

            /** \brief Sets the collision check caching to save calls to the collision
            checker with slightly memory usage as a counterpart */
            void setCacheCC(bool ccc)
            {
                cacheCC_ = ccc;
            }

            /** \brief Get the state of the collision check caching */
            bool getCacheCC() const
            {
                return cacheCC_;
            }

            void setExtended(bool efmt)
            {
                extendedFMT_ = efmt;
            }

            bool getExtended() const
            {
                return extendedFMT_;
            }

            void setOneSample(bool os)
            {
                oneSample_ = os;
            }

            bool getOneSample() const
            {
                return oneSample_;
            }

            void saveTree(const std::string &filename);

            // Specialized class for bi-directional trees
            class BiDirMotion {
            public:
                enum SetType { SET_NULL, SET_H, SET_W };
                
                BiDirMotion(TreeType* tree)
                    : state_(NULL), tree_(tree)
                {
                    parent_[FWD]        = NULL;
                    parent_[REV]        = NULL;
                    cost_[FWD]          = base::Cost(0.0);
                    cost_[REV]          = base::Cost(0.0);
                    hcost_[FWD]         = base::Cost(0.0);
                    hcost_[REV]         = base::Cost(0.0);
                    currentSet_[FWD]    = SET_W;
                    currentSet_[REV]    = SET_W;
                }

                /** \brief Constructor that allocates memory for the state */
                BiDirMotion(const base::SpaceInformationPtr &si, TreeType* tree)
                    : state_(si->allocState()), tree_(tree)
                {
                    parent_[FWD]        = NULL;
                    parent_[REV]        = NULL;
                    cost_[FWD]          = base::Cost(0.0);
                    cost_[REV]          = base::Cost(0.0);
                    hcost_[FWD]         = base::Cost(0.0);
                    hcost_[REV]         = base::Cost(0.0);
                    currentSet_[FWD]    = SET_W;
                    currentSet_[REV]    = SET_W;
                }

                typedef std::vector<BiDirMotion*> BiDirMotionPtrs;
                
                /** \brief The state contained by the motion */
                base::State             *state_;
                
                BiDirMotion*            parent_[2];         /**< The parent motion in the exploration tree */
                BiDirMotionPtrs         children_[2];       /**< The set of motions descending from the current motion */
                SetType                 currentSet_[2];
                TreeType*               tree_;
                base::Cost              cost_[2];           /**< The cost of this motion */

                /** \brief The minimum cost to go of this motion (heuristically computed) */
                base::Cost hcost_[2];
                /** \brief Contains the connections attempted FROM this node */
                std::set<BiDirMotion *> collChecksDone_;

                /** \brief Returns true if the connection to m has been already
                    tested and failed because of a collision */
                bool alreadyCC(BiDirMotion *m)
                {
                    if (collChecksDone_.find(m) == collChecksDone_.end())
                        return false;
                    return true;
                }

                /** \brief Caches a failed collision check to m */
                void addCC(BiDirMotion *m)
                {
                    collChecksDone_.insert(m);
                }

                /** \brief Set the cost to go heuristic cost */
                void setHeuristicCost(const base::Cost h)
                {
                    hcost_[*tree_] = h;
                }

                /** \brief Get the cost to go heuristic cost */
                base::Cost getHeuristicCost() const
                {
                    return hcost_[*tree_];
                }
                /** \brief Set the cost to go heuristic cost */
                void setOtherHeuristicCost(const base::Cost h)
                {
                    hcost_[(*tree_+1) % 2] = h;
                }

                /** \brief Get the cost to go heuristic cost */
                base::Cost getOtherHeuristicCost() const
                {
                    return hcost_[(*tree_+1) % 2];
                }

                inline base::Cost       getCost(void)           const           { return this->cost_[*tree_]; }
                inline base::Cost       getOtherCost(void)      const           { return this->cost_[(*tree_+1) % 2]; }
                inline void             setCost(double cost)                    { this->cost_[*tree_] = base::Cost(cost); }
                inline void             setCost(base::Cost cost)                { this->cost_[*tree_] = cost; }
                inline void             setOtherCost(base::Cost cost)            { this->cost_[(*tree_+1) % 2] = cost; }

                inline void             setParent(BiDirMotion* parent)          { this->parent_[*tree_] = parent; }
                inline void             setOtherParent(BiDirMotion* parent)     { this->parent_[(*tree_+1) % 2] = parent; }
                inline BiDirMotion*     getParent(void)         const           { return this->parent_[*tree_]; }
                inline BiDirMotion*     getAnyParent(int& tree)                 {
                    if (this->parent_[*tree_]) {
                        tree = 0;
                        return this->parent_[*tree_];
                    }
                    else if (this->parent_[(*tree_+1) % 2]) {
                        tree = 1;
                        return this->parent_[(*tree_+1) % 2];
                    }
                    else return NULL;
                }

                inline void             setChildren(BiDirMotionPtrs children)   { this->children_[*tree_] = children; }
                //TODO make these return &
                inline BiDirMotionPtrs  getChildren(void)       const           { return this->children_[*tree_]; }
                inline BiDirMotionPtrs  getOtherChildren(void)       const      { return this->children_[(*tree_+1) % 2]; }

                inline void             setCurrentSet(SetType set)              { this->currentSet_[*tree_] = set; }
                inline void             setOtherCurrentSet(SetType set)         { this->currentSet_[(*tree_+1) % 2] = set; }
                inline SetType          getCurrentSet(void)     const           { return this->currentSet_[*tree_]; }
                inline SetType          getOtherSet(void)       const           { return this->currentSet_[(*tree_+1) % 2]; }

                inline void             setTreeType(TreeType* treePtr)          { this->tree_       = treePtr; }
                inline TreeType         getTreeType(void)       const           { return *tree_; }

                /** \brief Set the state associated with the motion */
                void setState(base::State *state) {
                    state_ = state;
                }

                /** \brief Get the state associated with the motion */
                base::State* getState() const {
                    return state_;
                }
            };

            typedef std::vector<BiDirMotion*> BiDirMotionPtrs;

            struct BiDirMotionCompare {
                bool operator()(const BiDirMotion* p1, const BiDirMotion* p2) const {
                    if (heuristics_)
                        return opt_->isCostBetterThan(opt_->combineCosts(p1->getCost(), p1->getHeuristicCost()),
                                                      opt_->combineCosts(p2->getCost(), p2->getHeuristicCost()));
                    else
                        return opt_->isCostBetterThan(p1->getCost(), p2->getCost());
                }

                base::OptimizationObjective* opt_;
                bool heuristics_;
            };

            // Heap structure for storing the best node for expansion in either the fwd or rev trees
            typedef ompl::BinaryHeap<BiDirMotion*, BiDirMotionCompare> BiDirMotionBinHeap;
            
            // Tree functions
            void swapTrees();
            void useFwdTree()   { tree_ = FWD; }
            void useRevTree()   { tree_ = REV; }
            
            // Planning functions
            double distanceFunction(const BiDirMotion* a, const BiDirMotion* b) const {
                return opt_->motionCost(a->getState(), b->getState()).value();
            }
            
            double calculateUnitBallVolume(const unsigned int dimension) const;
            double calculateRadius(unsigned int dimension, unsigned int n) const;
            void freeMemory();
            void saveNeighborhood( boost::shared_ptr< NearestNeighbors<BiDirMotion*> > nn, BiDirMotion* m);
            
            void sampleFree( boost::shared_ptr<NearestNeighbors<BiDirMotion*> > nn,
                const base::PlannerTerminationCondition &ptc );

            void initializeProblem( base::GoalSampleableRegion*& goal_s );

            void expandTreeFromNode( BiDirMotion *&z, BiDirMotion *&connection_Point );
            
            bool plan( BiDirMotion *x_init, BiDirMotion *x_goal, BiDirMotion *&z, const base::PlannerTerminationCondition& ptc );

            // Set the termination condition
            bool termination( BiDirMotion *&z, BiDirMotion *&connection_Point, const base::PlannerTerminationCondition& ptc );
            
            // Choose exploration tree and node z to expand
            void chooseTreeAndExpansionNode( BiDirMotion *&z );
            
            // Trace a path of nodes along a tree towards the root (forward or reverse)
            void tracePath( BiDirMotion *z, BiDirMotionPtrs& path );

            void updateNeighborhood(BiDirMotion *m, const std::vector<BiDirMotion *> nbh);

            void updateKNeighborhood(BiDirMotion *m, const std::vector<BiDirMotion*> nbh);

            void insertNewSampleInOpen(const base::PlannerTerminationCondition& ptc);

            void steerNewSampleInOpen(const base::PlannerTerminationCondition& ptc);

            // Member variables
            unsigned int                numSamples_;
            double                      radiusMultiplier_;
            double                      freeSpaceVolume_;
            unsigned int                collisionChecks_;
            bool                        nearestK_;
            double                      NNr;
            unsigned int                NNk;
            TreeType                    tree_;
            ExploreType                 exploration_;
            TerminateType               termination_;
            const bool                  precomputeNN_;
            boost::shared_ptr< NearestNeighbors<BiDirMotion*> >     nn_;
            std::map<BiDirMotion*, BiDirMotionPtrs >                neighborhoods_;
            BiDirMotionBinHeap                                      H[2];           /**< Frontier */
            std::map<BiDirMotion*, BiDirMotionBinHeap::Element*>    H_elements[2];  /**< Frontier elements */

            /** \brief State sampler */
            base::StateSamplerPtr sampler_;

            /** \brief The cost objective function */
            base::OptimizationObjectivePtr opt_;

           /** \brief Flag to activate the cost to go heuristics */
            bool heuristics_;
            /** \brief Goal state caching to accelerate cost to go heuristic computation */
            base::State* heurGoalState_[2];

            bool cacheCC_;

            bool extendedFMT_;

            bool oneSample_;

            // For sorting a list of costs and getting only their sorted indices
            struct CostIndexCompare
            {
                CostIndexCompare(const std::vector<base::Cost>& costs,
                                 const base::OptimizationObjective &opt) :
                    costs_(costs), opt_(opt)
                {}
                bool operator()(unsigned i, unsigned j)
                {
                    return opt_.isCostBetterThan(costs_[i],costs_[j]);
                }
                const std::vector<base::Cost>& costs_;
                const base::OptimizationObjective &opt_;
            };

        };
        
    }   // End "geometric" namespace
}       // End "ompl" namespace


#endif	/* BIDIRECTIONALFMT_H */

