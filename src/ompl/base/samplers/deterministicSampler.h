#ifndef __FMT__DETERMINISTIC_STATE_SAMPLER_H__
#define __FMT__DETERMINISTIC_STATE_SAMPLER_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/StateSampler.h>
#include <omplapp/apps/SE3RigidBodyPlanning.h>

namespace ompl
{
    namespace base
    {
		
		 /// @cond IGNORE
        OMPL_CLASS_FORWARD(StateSpace);
        /// @endcond

        /// @cond IGNORE
        /** \brief Forward declaration of ompl::base::StateSampler */
        OMPL_CLASS_FORWARD(DeterministicStateSampler);
        /// @endcond
        
        /** \brief Definition of a real vector state space sampler that allows
        to sample predefined states. */
        class DeterministicStateSampler :
                public ompl::base::StateSampler
        {

        protected:
            unsigned int samplingCount;
            std::vector<State*> samples;
            std::string filename;

        public:
            DeterministicStateSampler(const StateSpace *ss, unsigned int maxSamples)
                : StateSampler(ss)
                , samplingCount(0)
            {
                //createNewSamples(maxSamples);
				std::string filename = "/home/jvgomez/Downloads/haltonSamples.txt";
                loadSamplesFromFile(filename);
            }

            // create deterministic state sampler with existing file
            DeterministicStateSampler(const StateSpace *ss, const std::string& fileName)
                : StateSampler(ss)
                , samplingCount(0), filename(fileName)
            {
                loadSamplesFromFile(fileName);
            }
            
            void createNewSamples(unsigned int maxSamples)
            {
                /*StateSamplerPtr sampler = space_->allocDefaultStateSampler();
                for(unsigned int i = 0; i < maxSamples; i++)
                {
                    State* rstate = space_->allocState();
                    sampler->sampleUniform(rstate);
                    samples.push_back(rstate);
                }*/
                std::string filename = "/home/jvgomez/Downloads/haltonSamples.txt";
                loadSamplesFromFile(filename);
            }

            void loadSamplesFromFile(const std::string& fileName)
            {
                std::string line;
                std::ifstream sampleFile(fileName.c_str());
                std::cout << "Loading from file " << filename.c_str() << std::endl;
                assert(sampleFile.is_open());

                // read dimension file
                getline (sampleFile,line);
                unsigned int dimension = (unsigned int)atoi(line.c_str());

                // check if dimension of data in file is sufficient
                if(dimension >= this->space_->getDimension())
                {
                    while ( getline (sampleFile,line) )
                    {
                        std::istringstream in(line);     //make a stream for the line itself

                        // get

                        // create new state and add to samples
                        State* state =  space_->allocState();
//                        RealVectorStateSpace::StateType *rstate = static_cast<RealVectorStateSpace::StateType*>(state);
//                        for (unsigned int i = 0 ; i < dimension ; ++i)
//                        {
//                            double val;
//                            in >> val;
//                            rstate->values[i] = val;
//                        }
                        std::vector<double> reals;

                        // TODO: needs to be dimension + 1 due to quaternions
                        unsigned dimTemp = dimension;
                        if(dimension == 6) {
                            dimTemp = dimension + 1;
                        }
                        for (unsigned i = 0 ; i < dimTemp ; ++i)
                        {
                            double val;
                            in >> val;
                            // std::cout << val << ' ';
                            reals.push_back(val);
                        }
                        // std::cout << std::endl;
                        space_->copyFromReals(state, reals);
                        samples.push_back(state);
                    }
                    sampleFile.close();
                }
            }


            ~DeterministicStateSampler()
            {
                for(unsigned int i = 0; i < samples.size(); i++)
                {
                    this->space_->freeState(samples[i]);
                }
            }


            void saveSamplesToFile(const std::string& fileName)
            {
                filename = fileName;
                std::ofstream sampleFile(fileName.c_str());
                assert(sampleFile.is_open());

                //const unsigned int dimension = space_->getDimension();
                //sampleFile << dimension << std::endl;
                // sampleFile << "200.00 -100.00 70.57 0.00 0.00 0.00 1.00" << std::endl;

                for(unsigned int i =0; i < samples.size(); i++)
                {
                    std::vector<double> reals;
                    space_->copyToReals(reals, samples[i]);
                    for(unsigned int j = 0; j < reals.size(); j++)
                    {
                        sampleFile << reals[j] <<" ";
                    }
                    sampleFile << std::endl;
                }

                sampleFile.close();
            }



            /* Reset the state sampler. */
            void reset()
            {
                samplingCount = 0;
            }

            /* override sample Uniform method that is called by the planning algorithms */
            void sampleUniform(State *state)
            {
                if (samplingCount == samples.size())
                {
					// TODO: What to do in this case?
					assert(0);
                    createNewSamples(samples.size());
                    saveSamplesToFile(filename);
                    OMPL_INFORM("Created additional samples and augmented the sample file");
                }
                assert(samplingCount < samples.size());
                
                space_->copyState(state, samples[samplingCount]);
                
//                const unsigned int dim = space_->getDimension();                
//                RealVectorStateSpace::StateType *rstate = static_cast<RealVectorStateSpace::StateType*>(state);
//                for (unsigned int i = 0 ; i < dim ; ++i)
//                    rstate->values[i] = static_cast<RealVectorStateSpace::StateType*>(samples[samplingCount])->values[i];

                samplingCount++;
            }

            void sampleUniformNear(State *state, const State *near, const double distance)
            {
                assert(0); //not implemented
/*                assert(samplingCount < samples.size());
                const unsigned int dim = space_->getDimension();

                std::vector<double> rreals;
                std::vector<double> nreals;
                space_->copyToReals(nreals, near);
                
                for (unsigned int i = 0 ; i < dim ; ++i)
                    rreals[i] = rng_.uniformReal(nreals[i] - distance, nreals[i] + distance);
                
                space_->copyFromReals(state, rreals);
                space_->enforceBounds(state);
                */
            }

            void sampleGaussian(State *state, const State *mean, const double stdDev)
            {
                assert(0); //not implemented
            }
            
            double getInformedMeasure(const Cost &currentCost) const {
				return space_->getMeasure();
			}
			
			double getInformedMeasure(const Cost &minCost, const Cost &maxCost) const {
				return space_->getMeasure();
			}
			
			bool hasInformedMeasure() const {
				return false;
			}
        };
    }
}

#endif // __FMT__DETERMINNISTIC_STATE_SAMPLER_H__
