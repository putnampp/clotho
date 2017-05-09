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
#ifndef CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_BATCHED_HPP_
#define CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_BATCHED_HPP_

#ifdef DEBUG_MODE
#define DEBUGGING 0
#endif  // DEBUG_MODE

#include "qtlsim_logger.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "clotho/genetics/population_growth_toolkit.hpp"

#include "clotho/data_spaces/allele_space/allele_space_vector.hpp"
#include "clotho/data_spaces/allele_space/allele_generators.hpp"

#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_space_generator.hpp"

#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"
#include "clotho/data_spaces/free_space/free_space_accumulator.hpp"

#include "clotho/data_spaces/population_space/population_spaces.hpp"
#include "clotho/data_spaces/selection/selection.hpp"

#include "clotho/data_spaces/fitness/general_fitness.hpp"

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/bit_helper.hpp"

#include "clotho/data_spaces/task/thread_count_parameter.hpp"
#include "clotho/data_spaces/task/thread_pool2.hpp"

#include "clotho/data_spaces/offspring_generator/offspring_generators.hpp"

#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

#include "clotho/data_spaces/mutation/mutation_allocator.hpp"

#include <vector>

struct parallel_pipeline_batched {};

#ifdef USE_OFF_GEN_PIPELINE
#define OFF_GEN_METHOD clotho::genetics::off_gen_pipelined
#else
#define OFF_GEN_METHOD clotho::genetics::off_gen_default
#endif  // USE_OFF_GEN_PIPELINE

template < class RNG, class RealType, class BlockType, class SizeType >
class Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline_batched > {
public:
    typedef Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline_batched >             self_type;

    typedef RealType                    position_type;
    typedef RealType                    weight_type;
    typedef weight_type *               phenotype_type;
    typedef BlockType                   block_type;

    typedef RNG                         random_engine_type;

    typedef SizeType                    size_type;

    typedef clotho::genetics::thread_pool2< RNG >                                                   thread_pool_type;
    typedef clotho::genetics::AlleleSpace< position_type, size_type >                   allele_type;

#ifdef USE_VECTOR_ALIGNMENT
    typedef clotho::genetics::index_vector_alignment< unsigned int, block_type >    alignment_type;
#else
    typedef clotho::genetics::row_block_alignment< block_type >                     alignment_type;
#endif

    typedef clotho::genetics::trait_space_vector< weight_type >                     trait_space_type;
    typedef clotho::genetics::population_space< alignment_type, trait_space_type >  sequence_space_type;

    typedef clotho::genetics::free_space_accumulator_mt< block_type, size_type >                    free_space_type;
    typedef typename free_space_type::buffer_type                                                   free_buffer_type;

    typedef clotho::genetics::mutation_allocator< random_engine_type, size_type >                 mutation_alloc_type;

    typedef clotho::genetics::AlleleGenerator< random_engine_type, allele_type >                allele_generator_type;
    typedef clotho::genetics::TraitSpaceGenerator2< trait_space_type >                          trait_generator_type;

    typedef clotho::genetics::GeneralFitness                                                      fitness_type;

    typedef std::shared_ptr< ipopulation_growth_generator >                                         population_growth_generator_type;
    typedef std::shared_ptr< ipopulation_growth >                                                   population_growth_type;

    typedef clotho::genetics::batch_offspring_generator< random_engine_type, sequence_space_type, allele_type, fitness_type, trait_space_type, free_buffer_type, OFF_GEN_METHOD >           offspring_generator_type;
    typedef typename offspring_generator_type::parameter_type parameter_type;

    typedef typename offspring_generator_type::mutation_pool_type                               mutation_pool_type;
    typedef typename offspring_generator_type::mutation_distribution_type                       mutation_distribution_type;

    typedef typename offspring_generator_type::time_vector  time_vector;

    typedef std::map< std::string, time_vector> time_map;

    friend struct clotho::utility::state_getter< self_type >;

    Engine( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_trait_space( config )
        , m_fixed_traits( config )
        , m_generation( 0 )
        , m_pop_growth()
        , m_thread_count( config )
        , m_free_space( m_thread_count.m_tc + 1 )
        , m_worker_rng( NULL )
        , m_mut_alloc( rng, config )
        , trait_gen( config )
        , allele_gen( rng, config )
        , m_recomb_rate( config )
        , m_bias_rate( config )
        , m_main_off_gen( rng, &m_allele_space, &m_trait_space, config )
        , m_worker_off_gens( NULL )
    {
        population_growth_generator_type tmp  = population_growth_toolkit::getInstance()->get_tool( config );
        if( tmp ) {
            m_pop_growth = tmp->generate();
            if( m_pop_growth ) {
                m_pop_growth->log( std::cerr );
                std::cerr << std::endl;
            }
        } else {
            population_growth_toolkit::getInstance()->tool_configurations( config );
        }

        init(0, config);
    }

    size_t getGeneration() const {
        return m_generation;
    }

    void init( size_t aN, boost::property_tree::ptree & config ) {
        size_t pN = 0;
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        boost::property_tree::ptree lconfig;

        lconfig = config.get_child( "generations", lconfig);

        unsigned int intermed = lconfig.get< unsigned int>( "intermediates", 0 );

        lconfig.put("intermediates", intermed );
        config.put_child( "generations", lconfig );

        for( unsigned int i = 0; i < intermed + 2; ++i ) {
            sequence_space_type * ss = new sequence_space_type();
            m_pop.push_back( ss );

            fitness_type * f = new fitness_type( config );

            assert( f->isEmpty() );

            m_fits.push_back( f );

            mutation_pool_type * p = new mutation_pool_type();
            m_mut_pools.push_back( p );

            mutation_distribution_type * d = new mutation_distribution_type();
            m_mut_dists.push_back( d );
        }


        if( m_thread_count.m_tc > 0 ) {
            m_worker_rng = new random_engine_type * [ m_thread_count.m_tc ];
            m_worker_off_gens = new offspring_generator_type * [ m_thread_count.m_tc ];

            for( int offset = 0; offset < m_thread_count.m_tc; ++offset ) {
                m_worker_rng[ offset ] = new random_engine_type( (*m_rand)());
                m_worker_off_gens[ offset ] = new offspring_generator_type( m_worker_rng[ offset ], &m_allele_space, &m_trait_space, config );
            }
        }

#ifdef USE_ROW_VECTOR
        m_pop1.getSequenceSpace().fill_empty();
        m_pop1.getSequenceSpace().finalize();
#endif // USE_ROW_VECTOR
        ++m_generation;
    }

    void simulate( ) {
        // use the last population as the first population for the next round
        std::swap( m_pop[ 0 ], m_pop[ m_pop.size() - 1 ] );
        std::swap( m_fits[ 0 ], m_fits[ m_fits.size() - 1] );

        size_t all_size = buildProjectedPopulations( );

#ifdef DEBUGGING 
        BOOST_LOG_TRIVIAL( debug ) << "Generation " << m_generation << ": " << pN << " individuals; " << pM << " new alleles";
        BOOST_LOG_TRIVIAL( debug ) << "Free space: " << free_count << "; alleles: " << m_allele_space.size();
        BOOST_LOG_TRIVIAL( debug ) << "Rescaling child population to be: " << pN << " individuals x " << all_size << " alleles";
#endif // DEBUGGING

        // grow the free space buffer to a size relative to the last projected population
        m_free_space.resetBuffers( m_pop.back()->getMaxBlocks() );

//        bool is_all_neutral = m_allele_space.isAllNeutral();
        bool is_all_neutral = m_trait_space.isAllNeutral();

        thread_pool_type tpool( m_thread_count.m_tc );
        
        const size_t TC = m_thread_count.m_tc + 1;

        std::vector< parameter_type * > params;

        for( int j = 0; j < m_thread_count.m_tc; ++j ) {
            free_buffer_type tbuf = m_free_space.getThreadBuffer( j );
            m_worker_off_gens[ j ]->reset( tbuf );

            unsigned int prev_start = 0, prev_end = m_pop[ 0 ]->getIndividualCount();

            for( unsigned int p = 1; p < m_pop.size(); ++p ) {
                unsigned int pN = m_pop[ p ]->getIndividualCount();
                size_t BATCH_SIZE = (pN / TC) + ((pN % TC > 0) ? 1 : 0);

                unsigned int off_idx = j * BATCH_SIZE;
                parameter_type * param = new parameter_type( m_pop[ p - 1 ], m_pop[ p ], m_fits[ p - 1], m_fits[ p ], m_mut_pools[ p - 1], m_mut_dists[ p - 1], prev_start, prev_end, off_idx, off_idx + BATCH_SIZE, is_all_neutral );

                m_worker_off_gens[ j ]->addPopulation( param );
                
                params.push_back( param );

                prev_start = off_idx;
                prev_end = off_idx + BATCH_SIZE;
            }
        }

        m_main_off_gen.reset( m_free_space.getThreadBuffer( m_thread_count.m_tc ) );

        unsigned int prev_start = 0, prev_end = m_pop[ 0 ]->getIndividualCount();
        for( unsigned int p = 1; p < m_pop.size(); ++p ) {
            unsigned int pN = m_pop[ p ]->getIndividualCount();
            size_t BATCH_SIZE = (pN / TC) + ((pN % TC > 0) ? 1 : 0);
            unsigned int off_idx = m_thread_count.m_tc * BATCH_SIZE;

            parameter_type * param = new parameter_type( m_pop[ p - 1 ], m_pop[ p ], m_fits[ p - 1], m_fits[ p ], m_mut_pools[ p - 1], m_mut_dists[ p - 1], prev_start, prev_end, off_idx, pN, is_all_neutral );
            m_main_off_gen.addPopulation( param );
            params.push_back( param );

            prev_start = off_idx;
            prev_end = pN;
        }

        tpool.post_list( m_worker_off_gens, m_thread_count.m_tc );

        m_main_off_gen();

        tpool.sync();

        m_free_space.finalize(all_size);

        recordAlleles();

        while( ! params.empty() ) {
            parameter_type * p = params.back();
            params.pop_back();
            delete p;
        }

        ++m_generation;
    }

    sequence_space_type * getChildPopulation() const {
        return m_pop.back();
    }

    sequence_space_type * getParentPopulation() const {
        return m_pop[ 0 ];
    }

    void recordAlleles() {
        for( unsigned int i = 1; i < m_pop.size(); ++i ) {
            clotho::utility::add_value_array( free_sizes, m_free_space.free_size() );
            clotho::utility::add_value_array( var_sizes, m_free_space.variable_count() );
            //clotho::utility::add_value_array( fixed_sizes, m_free_space.fixed_size() );
            clotho::utility::add_value_array( fixed_sizes, m_fixed.size() );
        }
    }

    void buildTimeVectorLog( boost::property_tree::ptree & log, time_map & times ) {
        for( typename time_map::iterator it = times.begin(); it != times.end(); it++ ) {
            boost::property_tree::ptree starts, stops;

            for( typename time_vector::iterator iter = it->second.begin(); iter != it->second.end(); ++iter ) {
                clotho::utility::add_value_array( starts, iter->first );
                clotho::utility::add_value_array( stops, iter->second );
            }

            boost::property_tree::ptree tmp;
            tmp.put_child( "start", starts );
            tmp.put_child( "stop", stops );

            log.put_child( it->first, tmp );
        }
    }

    void getPerformanceResults( boost::property_tree::ptree & log ) {

        for( int i = 0; i < m_thread_count.m_tc; ++i ) {
            std::ostringstream oss;
            oss << "W_" << (i + 1);

            boost::property_tree::ptree mut, xo, ph, fx, ft;
            m_worker_off_gens[ i ]->record( xo, mut, ph, fx, ft );
            
            log.put_child( "performance.mutate." + oss.str(), mut );
            log.put_child( "performance.crossover." + oss.str(), xo );
            log.put_child( "performance.fixed." + oss.str(), fx );
            log.put_child( "performance.phenotype." + oss.str(), ph );
            log.put_child( "performance.fitness." + oss.str(), ft );
        }

        boost::property_tree::ptree mut, xo, ph, fx, ft;
        m_main_off_gen.record( xo, mut, ph, fx, ft );
        
        log.put_child( "performance.mutate.main", mut );
        log.put_child( "performance.crossover.main", xo );
        log.put_child( "performance.fixed.main", fx );
        log.put_child( "performance.phenotype.main", ph );
        log.put_child( "performance.fitness.main", ft );

        log.put_child( "memory.free_count", free_sizes );
        log.put_child( "memory.variable_count", var_sizes );
        log.put_child( "memory.fixed_count", fixed_sizes );
    }

    allele_type * getAlleleSpace() {
        return &m_allele_space;
    }

    virtual ~Engine() { 
        if( m_worker_rng != NULL ) {
            for( int i = 0; i < m_thread_count.m_tc; ++i ) {
                delete m_worker_rng[ i ];
            }
            delete [] m_worker_rng;
            m_worker_rng = NULL;
        }

        if( m_worker_off_gens != NULL ) {
            for( int i = 0; i < m_thread_count.m_tc; ++i ) {
                delete m_worker_off_gens[ i ];
            }

            delete [] m_worker_off_gens;
            m_worker_off_gens = NULL;
        }

        while( !m_pop.empty() ) {
            sequence_space_type * p = m_pop.back();
            m_pop.pop_back();
            delete p;

            mutation_pool_type * pool = m_mut_pools.back();
            m_mut_pools.pop_back();
            delete pool;

            mutation_distribution_type * dist = m_mut_dists.back();
            m_mut_dists.pop_back();
            delete dist;

            fitness_type * f = m_fits.back();
            m_fits.pop_back();
            delete f; 
        }
    }

protected:

    size_t buildProjectedPopulations( ) {
        // remove fixed alleles
        size_t all_count = m_pop[0]->getMaxAlleles();
        size_t free_count = removeFixedAlleles( m_pop[0] );

        typename free_space_type::base_type::iterator it = m_free_space.free_begin(), end = m_free_space.free_end();

        unsigned int gen_offset = m_generation * (m_pop.size() - 1);

        for( unsigned int i = 1; i < m_pop.size(); ++i ) {
            all_count = grow_child_population( m_pop[ i - 1 ], m_pop[ i ], m_fits[ i ], m_mut_pools[ i - 1], m_mut_dists[ i - 1], all_count, free_count, it, end, gen_offset++ );

            free_count = (end - it);
        }

        return all_count;
    }

    size_t grow_child_population( sequence_space_type * parent, sequence_space_type * off, fitness_type * fit, mutation_pool_type * mut_pool, mutation_distribution_type * mut_dist, size_t all_count, size_t free_count, typename free_space_type::base_type::iterator & it, typename free_space_type::base_type::iterator & end, unsigned int gen ) {
        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = parent->getIndividualCount();
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        size_type pM = m_mut_alloc.allocate( 2 * pN );        // generate the number of new mutations

        size_t all_size = child_max_alleles( all_count, free_count, pM );   // rescale allele space for child population given free space from parent population and new allele count (pM)
        off->grow( pN, all_size, m_trait_space.trait_count() );               // grow the child population accordingly

        fit->resize( pN );
        generate_child_mutations( off, mut_pool, mut_dist, pM, it, end, gen );

        return all_size;
    }

    void generate_child_mutations( sequence_space_type * off, mutation_pool_type * mut_pool, mutation_distribution_type * mut_dist, unsigned int N, typename free_space_type::base_type::iterator & it, typename free_space_type::base_type::iterator & end, unsigned int gen ) {
//        std::cerr << "Child population size: " << m_child->haploid_genome_count() << std::endl;
//        std::cerr << "Mutation count: " << N << std::endl;
//        std::cerr << "Free space: " << m_free_space.free_size() << std::endl;
//
        resetMutationEvents( mut_pool, mut_dist, N, off->haploid_genome_count() + 1 );
//        m_allele_space.alignNeutralToPopulation( off->getMaxBlocks() );

        boost::random::uniform_int_distribution< unsigned int > seq_gen( 0, off->haploid_genome_count() - 1);

        unsigned int i = 0;
        while( i < N ) {
            typename free_space_type::size_type all_idx = m_allele_space.size();

            if( it != end ) {
                all_idx = *it++;
            } else {
                m_allele_space.grow();
                m_trait_space.grow();
            }

            unsigned int seq_idx = seq_gen( *m_rand );

            assert( all_idx < off->getMaxAlleles() );

            mut_pool->at( i ) = all_idx;
            mut_dist->at( seq_idx + 1 ) += 1;

            allele_gen( m_allele_space, all_idx, gen );
            trait_gen(*m_rand, m_trait_space, all_idx );

            ++i;
        }

        // scan right to produce m_mut_pool relative index ranges
        for( unsigned int i = 2; i < mut_dist->size(); ++i ) {
            mut_dist->at( i ) += mut_dist->at( i - 1 );
        }
        
    }

    void resetMutationEvents( mutation_pool_type * mut_pool, mutation_distribution_type * mut_dist, unsigned int N, unsigned int S ) {
        while( mut_pool->size() < N ) {
            mut_pool->push_back(0);
        }

        while( mut_dist->size() < S ) {
            mut_dist->push_back(0);
        }

        std::fill( mut_dist->begin(), mut_dist->end(), 0 );
    }

/**
 * estimate the maximum number of alleles in the child
 *
 * N_parent - number of alleles in the parent population
 * F_parent - number of free alleles in the parent population
 * M_child  - number of new alleles to be added the child population
 */
    size_t child_max_alleles( size_t N_parent, size_t F_parent, size_t M_child ) const {
#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(info) << "Parent alleles: " << N_parent << "; Free: " << F_parent << "; New Alleles: " << M_child;
        std::cerr << "Parent alleles: " << N_parent << "; Free: " << F_parent << "; New Alleles: " << M_child << std::endl;
#endif  // DEBUGGING

        if( F_parent >= M_child ) {
            // if there are more free alleles in the parent generation
            // than there are new alleles to be added to the child generation
            // then do not adjust scale of the allele space
            return N_parent;
        } else {
            return N_parent + (M_child - F_parent);
        }
    }

    size_type removeFixedAlleles( sequence_space_type * ss ) {
        typedef typename free_space_type::iterator  fixed_iterator;
        typedef typename trait_space_type::iterator trait_iterator;


        fixed_iterator  fix_it = m_free_space.fixed_begin();
        fixed_iterator  fix_end = m_free_space.fixed_end();

        while( fix_it != fix_end ) {
            size_type  fixed_index = *fix_it++;

            ss->remove_fixed_allele( fixed_index );

            m_fixed.append( m_allele_space, fixed_index );
            
            trait_iterator tstart = m_trait_space.begin( fixed_index ), tend = m_trait_space.end( fixed_index );
            m_fixed_traits.append( tstart, tend );
        }

#ifdef DEBUGGING
        typedef typename free_space_type::iterator free_iterator;
        free_iterator fr_it = m_free_space.free_begin();
        free_iterator fr_end = m_free_space.free_end();

        unsigned int j = 0;
        while( fr_it != fr_end ) {
            size_type i = *fr_it++;

            if( !ss->freeColumn( i ) ) {
                assert(false);
            }
            ++j;
        }

        assert( j == m_free_space.free_size() );
#endif // DEBUGGING

        return m_free_space.free_size();
    }

    random_engine_type  * m_rand;

    allele_type          m_allele_space, m_fixed;

    std::vector< sequence_space_type * > m_pop;

    trait_space_type        m_trait_space, m_fixed_traits;

    std::vector< fitness_type * > m_fits;

    size_t                  m_generation;

    population_growth_type  m_pop_growth;

    thread_count_parameter  m_thread_count;
    free_space_type         m_free_space;

    random_engine_type **   m_worker_rng;
    mutation_alloc_type     m_mut_alloc;
    trait_generator_type    trait_gen;
    allele_generator_type   allele_gen;

    recombination_rate_parameter< double > m_recomb_rate;
    sequence_bias_parameter< double > m_bias_rate;

    std::vector< mutation_pool_type * > m_mut_pools;
    std::vector< mutation_distribution_type * >  m_mut_dists;

    offspring_generator_type m_main_off_gen;
    offspring_generator_type ** m_worker_off_gens;

    time_map m_fixed_times, m_mutate_times, m_xover_times, m_pheno_times;
    boost::property_tree::ptree free_sizes, var_sizes, fixed_sizes, m_fitness_times;
};

namespace clotho {
namespace utility {

template < class RNG, class RealType, class BlockType, class SizeType >
struct state_getter< Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline_batched > > {
    typedef Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline_batched >           object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree tr;
        state_getter< typename object_type::trait_space_type > tr_logger;
        tr_logger( tr, obj.m_trait_space );

        boost::property_tree::ptree fr;
        state_getter< typename object_type::free_space_type > free_logger;
        free_logger( fr, obj.m_free_space );

        boost::property_tree::ptree fx, alls;
        state_getter< typename object_type::allele_type > all_logger;
        all_logger( fx, obj.m_fixed );
        all_logger( alls, obj.m_allele_space );

        boost::property_tree::ptree c_pop;
        state_getter< typename object_type::sequence_space_type > pop_logger;
        pop_logger( c_pop, *(obj.m_pop.back()) );

        boost::property_tree::ptree ft;
        clotho::utility::add_value_array( ft, obj.m_fits.back()->begin(), obj.m_fits.back()->end() );
        c_pop.put_child( "fitness", ft );

        s.put_child( "free_space", fr );
        s.put_child( "allele_space", alls );
        s.put_child( "trait_space", tr );
        s.put_child( "fixed_alleles", fx );

        s.put_child( "child", c_pop );
    }
};

}
}

#endif  // CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_BATCHED_HPP_
