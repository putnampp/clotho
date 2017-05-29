//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef ENGINE_CUDA_HOST_RANDOM_HPP_
#define ENGINE_CUDA_HOST_RANDOM_HPP_

#include <algorithm>
#include <boost/property_tree/ptree.hpp>

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

//#include "clotho/cuda/analysis/pairwise_difference.hpp"

#include "clotho/utility/timer.hpp"

#include "clotho/cuda/data_spaces/allele_space/host_allele_space.hpp"
#include "clotho/cuda/data_spaces/trait_space/host_trait_space.hpp"
#include "clotho/cuda/data_spaces/population_space/host_population_space.hpp"

#include "clotho/cuda/mutation/host_mutate_generator.hpp"
#include "clotho/cuda/selection/host_fit_selection_generator.hpp"
#include "clotho/cuda/crossover/host_crossover_generator.hpp"
#include "clotho/cuda/free_space/host_free_space.hpp"
#include "clotho/cuda/phenotype/host_phenotype_translator.hpp"

#include "clotho/cuda/fitness/host_fitness_translator.hpp"
#include "clotho/cuda/data_spaces/population_space/sample_population.hpp"

#include "clotho/cuda/analysis/host_allele_frequency.hpp"
#include "clotho/cuda/analysis/host_sequence_weight.hpp"
//#include "clotho/cuda/analysis/pairwise_difference.hpp"
//
#include "clotho/random/seed_parameter.hpp"

#include "engine_def.hpp"

struct cuda_host_random {};

template < class RNG, class RealType, class BlockType, class SizeType >
class Engine< RNG, RealType, BlockType, SizeType, cuda_host_random > : public clotho::utility::iStateObject {
public:

    typedef RNG                                                     random_engine_type;
    typedef RealType                                                real_type;
    typedef BlockType                                               block_type;

    typedef HostAlleleSpace< real_type >                            allele_space_type;
    typedef HostTraitSpace< real_type >                             trait_space_type;
    typedef HostPopulationSpace< real_type, unsigned int >          population_space_type;

    typedef HostMutateGenerator                                     mutation_generator_type;

    typedef HostSelectionGenerator                                  selection_generator_type;
    typedef HostCrossoverGenerator< typename allele_space_type::location_type > crossover_generator_type;

    typedef HostPhenotypeTranslator                                 phenotype_generator_type;

    typedef HostFreeSpace                                           free_space_type;

    typedef HostFitnessTranslator< typename population_space_type::fitness_type > fitness_translator_type;

    typedef clotho::utility::timer                                  timer_type;

    typedef typename population_space_type::sequence_space_type::bit_helper_type bit_helper_type;

    Engine( boost::property_tree::ptree & config ) :
         prev_pop( &hPop1 )
        , current_pop( &hPop0 )
        , traits( config )
        , mut_gen( config )
        , xover_gen( config )
        , sel_gen(config)
        , fit_trans( config )
        , m_generation(0)
    {
        parse_config( config );

        clotho::cuda::curand_state_pool::getInstance()->initialize( config );

        initialize();
    }

    void simulate( unsigned int N ) {
        std::cerr << "Generation: " << m_generation << std::endl;

        swap();
        unsigned int cur_seq_count = 2 * N;

        // use fixed space to update fixed alleles
        recordFixed( prev_pop );

        // use free space to determine new allele count for current population
        unsigned int allele_count = mut_gen.initialize( m_rng, prev_pop, cur_seq_count, alleles.getAlleleCount() );

        std::cerr << "Allele count: " << allele_count << std::endl;
        // resize the host alleles
        alleles.resize( allele_count );
        // resize the host traits
        traits.resize( alleles );
        // resize the current population accordingly
        current_pop->resize( alleles, traits, cur_seq_count );

        xover_gen.initialize( m_rng, current_pop );
        xover_gen.buildRecombinationMask( alleles, current_pop );

        sel_gen( m_rng, prev_pop, current_pop );

        xover_gen.performCrossover( prev_pop, current_pop, sel_gen );

        mut_gen.generate( m_rng, alleles, traits, m_generation++ );

        mut_gen( current_pop );

        // evaluate the free space of the current population after mutation
        free_space_type fr;
        fr( current_pop );

        phenotype_generator_type ph;
        ph( current_pop, traits, alleles.getAlleleCount() );

        fit_trans( current_pop );
        
        // in this model, mutation occurs after recombination
        // so new mutation locations do impact recombination process
        // and updating the device's allele location vector can be delayed
        // push locations for next generation to device
        alleles.updateDevice();

        cudaDeviceSynchronize();

/*
        size_t fsize, tsize;
        cudaError_t err = cudaMemGetInfo( &fsize, &tsize );

        if( err != cudaSuccess ) {
            std::cerr << "Unable to determine device memory (simulate): " << cudaGetErrorString(err) << std::endl;

            assert(false);
            fsize = 0;
        }
        clotho::utility::add_value_array(dev_space.get_child("free"), fsize );
*/
    }

    void trackStats() {
        clotho::utility::add_value_array( m_allele_track, alleles.getAlleleCount() );
        clotho::utility::add_value_array( m_free_track, prev_pop->getFreeCount() );
        clotho::utility::add_value_array( m_fixed_track, fixed_alleles.getAlleleCount());

        std::cerr << alleles.getAlleleCount() << ", " << prev_pop->getFreeCount() << ", " << fixed_alleles.getAlleleCount() << std::endl;
    }

    void recordFixed( population_space_type * pop ) {
        // update host downloads free/fixed space from device to host
        pop->updateHost();

        if( pop->getFixedCount() == 0 ) return;

        block_type * fixed = pop->getFreeSpace();
        unsigned int N = pop->getBlocksPerSequence();

        if( N > 0 ) {
            // trigger kernel to remove fixed alleles from population
            //
            free_space_type fr;
            fr.perform_remove_fixed( pop );
        }

        std::set< unsigned int > offsets;
        for( unsigned int i = 0; i < N; ++i ) {
            block_type b = fixed[ i ];

            if( b ) {
                for( unsigned int k = 0; k < bit_helper_type::BITS_PER_BLOCK; k++ ) {
                    if( b & 1) {
                        offsets.insert( i * bit_helper_type::BITS_PER_BLOCK + k );
                    }

                    b >>= 1;
                }
            }
        }

        if(!offsets.empty() ) {
            fixed_alleles.resize( fixed_alleles.getAlleleCount() + offsets.size() );

            for( std::set< unsigned int >::iterator it = offsets.begin(); it != offsets.end(); ++it ) {
                unsigned int idx = *it;

                fixed_alleles.push_back( alleles, idx );
            }
        }
    }

    population_space_type   * get_parent_population() {
        return prev_pop;
    }

    population_space_type   * get_offspring_population() {
        return current_pop;
    }

    void analyze_population( ) {
        all_freq( current_pop, alleles.getAlleleCount() );
        seq_weight( current_pop );
//
//        //pair_diff.evaluate( current_pop );
//
//        cudaDeviceSynchronize();
    }
//
    void analyze_sample( unsigned int N, boost::property_tree::ptree & samp_res) {
//        clotho::utility::timer t;
//
//        SamplePopulation< population_space_type > samps( current_pop, N );
//
////        pair_diff.evaluate( current_pop, samps.m_dSubpop );
//
//        cudaDeviceSynchronize();
//        t.stop();
//
////        pair_diff.get_state( samp_res );
//        samp_res.put( "rt", t.elapsed_long());
    }

    void get_tracked_states( boost::property_tree::ptree & state ) {
        state.put_child( "allele_space.sizes", m_allele_track );
        state.put_child( "free_space.sizes", m_free_track );
        state.put_child( "fixed_space.sizes", m_fixed_track );
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree cur, prev;
        current_pop->get_state( cur );
        prev_pop->get_state( prev );

        boost::property_tree::ptree al;
        alleles.get_state( al );

        boost::property_tree::ptree fx;
        fixed_alleles.get_state( fx );

        boost::property_tree::ptree sel;
        sel_gen.get_state( sel );

        boost::property_tree::ptree asis;
        all_freq.get_state( asis );
        seq_weight.get_state( asis );
////        pair_diff.get_state( asis );

//        state.put_child("device.memory", dev_space );

        state.put_child( "population.current", cur );
        state.put_child( "population.previous", prev );
        state.put_child( "selection", sel );
        state.put_child( "analysis", asis );
        state.put_child( "alleles.population", al );
        state.put_child( "alleles.fixed", fx );
    }

    void swap() {
        std::swap( prev_pop, current_pop );
    }

    virtual ~Engine() { }

protected:

    void parse_config( boost::property_tree::ptree & config ) {
        seed_parameter<> seed( config );

        m_rng.seed( seed.m_seed );
    }

    void initialize( ) {
/*
        size_t fsize = 0;
        size_t tsize = 0;

        cudaError_t err = cudaMemGetInfo( &fsize, &tsize);

        if( err != cudaSuccess ) {
            std::cerr << "Unable to initialize memory space sizes" << std::endl;
        }

        boost::property_tree::ptree tmp;
        clotho::utility::add_value_array(tmp, fsize );
        dev_space.put_child( "free", tmp );
        dev_space.put("total", tsize);
*/
    }

    random_engine_type      m_rng;
    population_space_type   hPop0, hPop1;
    population_space_type   * prev_pop, * current_pop;

    allele_space_type          alleles, fixed_alleles;

    trait_space_type            traits;

    mutation_generator_type     mut_gen;
    crossover_generator_type    xover_gen;
    selection_generator_type    sel_gen;

    fitness_translator_type     fit_trans;

    HostAlleleFrequency         all_freq;
    HostSequenceWeight          seq_weight;
//    PairwiseDifference          pair_diff;

    unsigned int                m_generation;
//    boost::property_tree::ptree dev_space;
//
    boost::property_tree::ptree m_allele_track, m_free_track, m_fixed_track;
};

#endif  // ENGINE_CUDA_HOST_RANDOM_HPP_
