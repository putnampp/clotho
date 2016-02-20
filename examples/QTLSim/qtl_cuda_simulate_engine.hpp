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
#ifndef QTL_CUDA_SIMULATE_ENGINE_HPP_
#define QTL_CUDA_SIMULATE_ENGINE_HPP_

#include <algorithm>
#include <boost/property_tree/ptree.hpp>

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"
#include "clotho/cuda/mutation/mutation_event_generator.hpp"

#ifndef USE_OFFSPRING_CROSSOVER
#include "clotho/cuda/crossover/crossover_generator.hpp"
#else
#include "clotho/cuda/crossover/offspring_generator.hpp"
#endif  // USE_OFFSPRING_CROSSOVER

//#include "clotho/cuda/selection/selection_event_generator.hpp"
#include "clotho/cuda/selection/fit_selection_generator.hpp"
#include "clotho/cuda/phenotype/phenotype_translator.hpp"

#include "clotho/cuda/analysis/allele_frequency.hpp"
#include "clotho/cuda/analysis/sequence_hamming_weight.hpp"
#include "clotho/cuda/analysis/pairwise_difference.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/cuda/fitness/quadratic_fitness_translator.hpp"
#include "clotho/cuda/data_spaces/population_space/sample_population.hpp"

template < class PopulationSpaceType > 
class qtl_cuda_simulate_engine : public clotho::utility::iStateObject {
public:
    typedef PopulationSpaceType                                 population_space_type;

    typedef typename population_space_type::allele_space_type       allele_space_type;
    typedef typename population_space_type::sequence_space_type     sequence_space_type;

    typedef typename allele_space_type::real_type                           real_type;
    typedef typename sequence_space_type::int_type                          int_type;
    typedef typename population_space_type::free_space_type::order_tag_type order_tag_type;

    typedef MutationEventGenerator< real_type, int_type, order_tag_type >   mutation_generator_type;
    typedef typename mutation_generator_type::space_type            mutation_event_space_type;

#ifndef USE_OFFSPRING_CROSSOVER
    typedef CrossoverGenerator< population_space_type >             crossover_generator_type;
#else
    typedef OffspringGenerator< population_space_type >             crossover_generator_type;
#endif  // USE_OFFSPRING_CROSSOVER

    //typedef SelectionEventGenerator                                 selection_generator_type;
    typedef FitSelectionGenerator< int_type, real_type >            selection_generator_type;

    typedef QuadraticFitnessTranslator< typename population_space_type::phenotype_space_type >                 fitness_type;

    typedef clotho::utility::timer                                  timer_type;

    qtl_cuda_simulate_engine( boost::property_tree::ptree & config ) :
        current_pop( &hPop0 )
        , prev_pop( &hPop1 )
        , mut_gen( config )
        , xover_gen( config )
        , sel_gen(config)
        , pheno_trans(config)
        , fit_trans( config )
        , all_freq( config )
        , seq_weight( config )
        , pair_diff( config )
    {
        parse_config( config );
        initialize();
    }

#ifndef USE_OFFSPRING_CROSSOVER
    void simulate( unsigned int N ) {
        swap();
        unsigned int cur_seq_count = 2 * N;

        prev_pop->record_fixed( fixed_alleles );

        mut_gen.generate( dMutations, cur_seq_count );

        current_pop->resize( prev_pop, dMutations, cur_seq_count );

        xover_gen( current_pop );

        // recombine parents 
        sel_gen.generate_and_recombine( prev_pop, current_pop, fit_trans.get_device_space() );

        mut_gen.scatter( current_pop, dMutations, cur_seq_count );
        
        pheno_trans.translate( current_pop );

        fit_trans( current_pop->pheno_space, cur_seq_count );
        
        current_pop->update_metadata();

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
#else 
    void simulate( unsigned int N ) {
        swap();
        unsigned int cur_seq_count = 2 * N;

        prev_pop->record_fixed( fixed_alleles );

        mut_gen.generate( dMutations, cur_seq_count );

        current_pop->resize( prev_pop, dMutations, cur_seq_count );

        sel_gen.generate( prev_pop, current_pop, fit_trans.get_device_space() );

        xover_gen( prev_pop, sel_gen.get_device_space(), current_pop );

        mut_gen.scatter( current_pop, dMutations, cur_seq_count );
        
        pheno_trans.translate( current_pop );

        fit_trans( current_pop->pheno_space, cur_seq_count );
        
        current_pop->update_metadata();

        cudaDeviceSynchronize();
    }
#endif  // USE_OFFSPRING_CROSSOVER

    void analyze_population( ) {
        all_freq.evaluate( current_pop );
        seq_weight.evaluate( current_pop );

        pair_diff.evaluate( current_pop );

        cudaDeviceSynchronize();
    }

    void analyze_sample( unsigned int N, boost::property_tree::ptree & samp_res) {
        clotho::utility::timer t;

        SamplePopulation< population_space_type > samps( current_pop, N );

        pair_diff.evaluate( current_pop, samps.m_dSubpop );

        cudaDeviceSynchronize();
        t.stop();

        pair_diff.get_state( samp_res );
        samp_res.put( "rt", t.elapsed_long());
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree cur, prev;
        current_pop->get_state( cur );
        prev_pop->get_state( prev );

        boost::property_tree::ptree fx;
        fixed_alleles.get_state( fx );

        boost::property_tree::ptree mut;
        get_device_object_state( mut, dMutations );

        boost::property_tree::ptree sel;
        sel_gen.get_state( sel );

        boost::property_tree::ptree asis;
        all_freq.get_state( asis );
        seq_weight.get_state( asis );
        pair_diff.get_state( asis );

        boost::property_tree::ptree fit;
        fit_trans.get_state( fit );

//        state.put_child("device.memory", dev_space );

        state.put_child( "population.current", cur );
        state.put_child( "population.current.fitness", fit );

        state.put_child( "population.previous", prev );
        state.put_child( "mutations", mut );
        state.put_child( "selection", sel );
        state.put_child( "analysis", asis );
        state.put_child( "alleles.fixed", fx );
    }

    void swap() {
        std::swap( prev_pop, current_pop );
    }

    virtual ~qtl_cuda_simulate_engine() {
        delete_space( dMutations );
    }

protected:

    void parse_config( boost::property_tree::ptree & config ) {
    }

    void initialize( ) {
        create_space( dMutations );

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

    mutation_event_space_type     * dMutations;
    population_space_type   hPop0, hPop1;
    population_space_type   * prev_pop, * current_pop;

    allele_space_type       fixed_alleles;

    mutation_generator_type     mut_gen;
    crossover_generator_type    xover_gen;
    selection_generator_type    sel_gen;
    phenotype_translator        pheno_trans;

    fitness_type                fit_trans;

    AlleleFrequency             all_freq;
    SequenceHammingWeight       seq_weight;
    PairwiseDifference          pair_diff;

//    boost::property_tree::ptree dev_space;
};

#endif  // QTL_CUDA_SIMULATE_ENGINE_HPP_
