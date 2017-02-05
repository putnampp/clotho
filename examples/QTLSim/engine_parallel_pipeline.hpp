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
#ifndef CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_HPP_
#define CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_HPP_

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

#include "clotho/data_spaces/task/thread_count_parameter.hpp"
#include "clotho/data_spaces/task/thread_pool2.hpp"

#include "clotho/data_spaces/offspring_generator/offspring_generators.hpp"

#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

#include "clotho/data_spaces/mutation/mutation_allocator.hpp"

#include <vector>

struct parallel_pipeline {};

template < class RNG, class RealType, class BlockType, class SizeType >
class Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline > {
public:
    typedef Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline >             self_type;

    typedef RealType                    position_type;
    typedef RealType                    weight_type;
    typedef weight_type *               phenotype_type;
    typedef BlockType                   block_type;

    typedef RNG                         random_engine_type;

    typedef SizeType                    size_type;

    typedef clotho::genetics::thread_pool2< RNG >                                                   thread_pool_type;
    typedef clotho::genetics::AlleleSpace< position_type, size_type >                               allele_type;

#ifdef USE_ROW_MODIFICATION
    typedef clotho::genetics::population_space_row_modified< block_type, weight_type >              sequence_space_type;
#else
    typedef clotho::genetics::population_space_row< block_type, weight_type >                       sequence_space_type;
#endif  // USE_ROW_MODIFICATION
    typedef clotho::genetics::trait_space_vector< weight_type >                                     trait_space_type;

    typedef clotho::genetics::free_space_accumulator_mt< block_type, size_type >                    free_space_type;
    typedef typename free_space_type::buffer_type                                                   free_buffer_type;

    typedef clotho::genetics::mutation_allocator< random_engine_type, size_type >                 mutation_alloc_type;

    typedef clotho::genetics::AlleleGenerator< random_engine_type, allele_type >                allele_generator_type;
    typedef clotho::genetics::TraitSpaceGenerator2< trait_space_type >                          trait_generator_type;


    typedef clotho::genetics::GeneralFitness                                                      fitness_type;
    typedef clotho::genetics::SelectionGenerator< random_engine_type, clotho::genetics::fitness_selection< fitness_type > >          selection_type;

    typedef std::shared_ptr< ipopulation_growth_generator >                                         population_growth_generator_type;
    typedef std::shared_ptr< ipopulation_growth >                                                   population_growth_type;

    typedef clotho::genetics::offspring_generator< random_engine_type, sequence_space_type, allele_type, selection_type, trait_space_type, free_buffer_type >           offspring_generator_type;

    typedef typename offspring_generator_type::mutation_pool_type                               mutation_pool_type;
    typedef typename offspring_generator_type::mutation_distribution_type                       mutation_distribution_type;
    typedef typename sequence_space_type::bit_helper_type bit_helper_type;

    friend struct clotho::utility::state_getter< self_type >;

    Engine( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_parent( &m_pop0 )
        , m_child( &m_pop1 )
        , m_trait_space( config )
        , m_fixed_traits( config )
        , select_gen( rng, config )
        , m_fit( config )
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

        init(0);
    }

    size_t getGeneration() const {
        return m_generation;
    }

    void init( size_t aN ) {
        size_t pN = 0;
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        m_pop1.grow( pN, aN, m_trait_space.trait_count() );
        m_pop1.clear();
        m_pop0.clear();

        if( m_thread_count.m_tc > 0 ) {
            m_worker_rng = new random_engine_type * [ m_thread_count.m_tc ];
            for( int i = 0; i < m_thread_count.m_tc; ++i ) {
                m_worker_rng[ i ] = new random_engine_type( (*m_rand)());
            }
        }

#ifdef USE_ROW_VECTOR
        m_pop1.getSequenceSpace().fill_empty();
        m_pop1.getSequenceSpace().finalize();
#endif // USE_ROW_VECTOR
        ++m_generation;
    }

    void simulate( ) {
        std::swap( m_child, m_parent );     // use the current child population as the parent population for the next round

        // remove fixed alleles
        size_t free_count = removeFixedAlleles( m_parent );

        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = select_gen.individual_count();
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        size_type pM = m_mut_alloc.allocate( 2 * pN );        // generate the number of new mutations

        size_t all_size = child_max_alleles( m_allele_space.size(), free_count, pM );   // rescale allele space for child population given free space from parent population and new allele count (pM)

#ifdef DEBUGGING 
        BOOST_LOG_TRIVIAL( debug ) << "Generation " << m_generation << ": " << pN << " individuals; " << pM << " new alleles";
        BOOST_LOG_TRIVIAL( debug ) << "Free space: " << free_count << "; alleles: " << m_allele_space.size();
        BOOST_LOG_TRIVIAL( debug ) << "Rescaling child population to be: " << pN << " individuals x " << all_size << " alleles";
#endif // DEBUGGING

        m_child->grow( pN, all_size, m_trait_space.trait_count() );               // grow the child population accordingly

        select_gen.update( m_fit, pN );

        m_mut_pool.clear();
        m_mut_dist.clear();

        generate_child_mutations( pM );
        thread_pool_type tpool( m_thread_count.m_tc );

        block_type * neutrals = new block_type[ m_child->getMaxBlocks() ];
        memset( neutrals, 0, m_child->getMaxBlocks() * sizeof( block_type ) );

        bool allNeutral = fill_neutrals( neutrals );

        const size_t TC = tpool.pool_size() + 1;
        const size_t BATCH_SIZE = (pN / TC) + ((pN % TC > 0)? 1:0);

        std::vector< offspring_generator_type * > task_list;
        task_list.reserve( TC );

        m_free_space.resetBuffers( m_child->getMaxBlocks() );

        unsigned int off_idx = 0, j = 0;
        while( off_idx + BATCH_SIZE < pN ) {
            free_buffer_type tbuf = m_free_space.getThreadBuffer( j );
            offspring_generator_type * ogen = new offspring_generator_type( m_worker_rng[ j ], m_parent, m_child, &m_allele_space, &m_mut_pool, &m_mut_dist, &select_gen, &m_trait_space, neutrals, tbuf, off_idx, off_idx + BATCH_SIZE, m_recomb_rate.m_rho, m_bias_rate.m_bias, allNeutral);
            task_list.push_back( ogen );
            off_idx += BATCH_SIZE;
            ++j;
        }

        tpool.post_list( task_list );

        if( off_idx < pN ) {
            free_buffer_type tbuf = m_free_space.getThreadBuffer( j );
            offspring_generator_type ogen( m_rand, m_parent, m_child, &m_allele_space, &m_mut_pool, &m_mut_dist, &select_gen, &m_trait_space, neutrals, tbuf, off_idx, pN, m_recomb_rate.m_rho, m_bias_rate.m_bias, allNeutral );
            ogen();

            // replace Allele space
            recordTimes( &ogen, "main" );
        }

        tpool.sync();

        while( !task_list.empty() ) {
            std::ostringstream oss;

            oss << "W_" << task_list.size();

            offspring_generator_type * tmp = task_list.back();
            task_list.pop_back();

            recordTimes( tmp, oss.str() );

            delete tmp;
        }

        timer_type fit_time;
        m_fit( m_child );
        fit_time.stop();

        clotho::utility::add_value_array( m_fitness_times, fit_time );

        m_free_space.finalize(all_size);

        recordAlleles();

        delete [] neutrals;

        ++m_generation;
    }

    bool fill_neutrals( block_type * neutrals ) {
        typename allele_type::neutral_vector::size_type i = m_allele_space.getNeutrals().find_first();

        while( i != allele_type::neutral_vector::npos ) {
            typename allele_type::neutral_vector::size_type idx = i / bit_helper_type::BITS_PER_BLOCK;
            typename allele_type::neutral_vector::size_type offset = i % bit_helper_type::BITS_PER_BLOCK;

            neutrals[ idx ] |= ((block_type) 1 << offset );
            
            i = m_allele_space.getNeutrals().find_next(i);
        }

        return m_allele_space.getNeutrals().all();
    }

    sequence_space_type * getChildPopulation() const {
        return m_child;
    }

    sequence_space_type * getParentPopulation() const {
        return m_parent;
    }

    void recordTimes( offspring_generator_type * gen, std::string key ) {

        boost::property_tree::ptree xover_times, mutate_times, pheno_times, fixed_times;

        xover_times = m_xover_times.get_child( key, xover_times );
        mutate_times = m_mutate_times.get_child( key, mutate_times );
        pheno_times = m_pheno_times.get_child( key, pheno_times );
        fixed_times = m_fixed_times.get_child( key, fixed_times );

        gen->recordCrossover( xover_times );
        gen->recordMutate( mutate_times );
        gen->recordPhenotype( pheno_times );
        gen->recordFixed( fixed_times );

        m_xover_times.put_child( key, xover_times );
        m_mutate_times.put_child( key, mutate_times );
        m_pheno_times.put_child( key, pheno_times );
        m_fixed_times.put_child( key, fixed_times );
    }

    void recordAlleles() {
        clotho::utility::add_value_array( free_sizes, m_free_space.free_size() );
        clotho::utility::add_value_array( var_sizes, m_free_space.variable_count() );
        clotho::utility::add_value_array( fixed_sizes, m_free_space.fixed_size() );
    }

    void getPerformanceResults( boost::property_tree::ptree & log ) {
        log.put_child( "performance.mutate", m_mutate_times );
        log.put_child( "performance.crossover", m_xover_times );
        log.put_child( "performance.fixed", m_fixed_times );
        log.put_child( "performance.phenotypes", m_pheno_times );
        log.put_child( "performance.fitness", m_fitness_times );

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
    }

protected:

    void generate_child_mutations( unsigned int N ) {
//        std::cerr << "Child population size: " << m_child->haploid_genome_count() << std::endl;
//        std::cerr << "Mutation count: " << N << std::endl;
//        std::cerr << "Free space: " << m_free_space.free_size() << std::endl;
//
        resetMutationEvents( N );

        boost::random::uniform_int_distribution< unsigned int > seq_gen( 0, m_child->haploid_genome_count() - 1);
        typename free_space_type::base_type::iterator it = m_free_space.free_begin(), end = m_free_space.free_end();

        unsigned int i = 0;
        while( i < N && it != end ) {
            typename free_space_type::size_type all_idx = *it++;
            unsigned int seq_idx = seq_gen( *m_rand );

            m_mut_pool[ i++ ] = all_idx;
            m_mut_dist[ seq_idx + 1 ] += 1;

            allele_gen( m_allele_space, all_idx, m_generation );
            trait_gen(*m_rand, m_trait_space, all_idx );
        }

        while( i < N ) {
            typename free_space_type::size_type all_idx = m_allele_space.size();
            unsigned int seq_idx = seq_gen( *m_rand );

            m_mut_pool[i++] = all_idx;
            m_mut_dist[ seq_idx + 1 ] += 1;

            assert( all_idx < m_child->getMaxAlleles() );

            allele_gen( m_allele_space, all_idx, m_generation );
            trait_gen( *m_rand, m_trait_space, all_idx );
            --N;
        }

        // scan right to produce m_mut_pool relative index ranges
        for( unsigned int i = 2; i < m_mut_dist.size(); ++i ) {
            m_mut_dist[ i ] += m_mut_dist[ i - 1 ];
        }
    }

    void resetMutationEvents( unsigned int N ) {
        m_mut_pool.clear();
        m_mut_pool.reserve( N );

        while( m_mut_pool.size() < N ) {
            m_mut_pool.push_back(0);
        }

        m_mut_dist.reserve( m_child->haploid_genome_count() + 1 );

        unsigned int i = 0;
        while( i < m_mut_dist.size() ) {
            m_mut_dist[ i ] = 0;
            ++i;
        }

        while( i < m_child->haploid_genome_count() + 1 ) {
            m_mut_dist.push_back(0);
            ++i;
        }
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

    allele_type             m_allele_space, m_fixed;
    sequence_space_type  m_pop0, m_pop1;
    sequence_space_type  * m_parent, * m_child;

    trait_space_type        m_trait_space, m_fixed_traits;


    selection_type          select_gen;
    fitness_type            m_fit;

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

    mutation_pool_type  m_mut_pool;
    mutation_distribution_type  m_mut_dist;
//
    boost::property_tree::ptree m_fixed_times, m_mutate_times, m_xover_times, m_pheno_times, m_fitness_times;
    boost::property_tree::ptree free_sizes, var_sizes, fixed_sizes;
};

namespace clotho {
namespace utility {

template < class RNG, class RealType, class BlockType, class SizeType >
struct state_getter< Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline > > {
    typedef Engine< RNG, RealType, BlockType, SizeType, parallel_pipeline >           object_type;

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
        pop_logger( c_pop, *(obj.m_child) );

        boost::property_tree::ptree ft;
        clotho::utility::add_value_array( ft, obj.m_fit.begin(), obj.m_fit.end() );
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

#endif  // CLOTHO_SIM_ENGINE_PARALLEL_PIPELINE_HPP_
