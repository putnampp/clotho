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
#include "clotho/utility/bit_helper.hpp"

#include "clotho/data_spaces/task/thread_count_parameter.hpp"
#include "clotho/data_spaces/task/thread_pool2.hpp"

#include "clotho/data_spaces/offspring_generator/offspring_generators.hpp"

#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

#include "clotho/data_spaces/mutation/mutation_allocator.hpp"

#include "clotho/initializer/initializers.hpp"

#include <vector>
#include <algorithm>

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

    typedef clotho::genetics::thread_pool2< RNG >                                   thread_pool_type;
    typedef clotho::genetics::AlleleSpace< position_type, size_type >               allele_type;

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

    typedef clotho::genetics::offspring_generator< random_engine_type, sequence_space_type, allele_type, fitness_type, trait_space_type, free_buffer_type >           offspring_generator_type;

    typedef typename offspring_generator_type::mutation_pool_type                               mutation_pool_type;
    typedef typename offspring_generator_type::mutation_distribution_type                       mutation_distribution_type;

    typedef typename offspring_generator_type::time_vector  time_vector;

    typedef clotho::genetics::PopulationInitializer population_init_type;

    typedef std::map< std::string, time_vector> time_map;

    friend struct clotho::utility::state_getter< self_type >;

    Engine( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_parent( &m_pop0 )
        , m_child( &m_pop1 )
        , m_trait_space( config )
        , m_fixed_traits( config )
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
        , m_main_off_gen( rng, &m_allele_space, &m_mut_pool, &m_mut_dist, &m_fit, &m_trait_space, m_recomb_rate.m_rho, m_bias_rate.m_bias )
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

        population_init_type pinit(config);

        init(pinit);
    }

    size_t getGeneration() const {
        return m_generation;
    }

    void init( population_init_type & pinit ) {
        size_t pN = 0;
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        if( m_thread_count.m_tc > 0 ) {
            m_worker_rng = new random_engine_type * [ m_thread_count.m_tc ];
            m_worker_off_gens = new offspring_generator_type * [ m_thread_count.m_tc ];
            for( int i = 0; i < m_thread_count.m_tc; ++i ) {
                m_worker_rng[ i ] = new random_engine_type( (*m_rand)());
                m_worker_off_gens[i] = new offspring_generator_type( m_worker_rng[ i ], &m_allele_space, &m_mut_pool, &m_mut_dist, &m_fit, &m_trait_space, m_recomb_rate.m_rho, m_bias_rate.m_bias );
            }
        }

#ifdef USE_ROW_VECTOR
        m_pop1.getSequenceSpace().fill_empty();
        m_pop1.getSequenceSpace().finalize();
#endif // USE_ROW_VECTOR

        buildInitialPopulation( pinit, pN );
        
        ++m_generation;
    }

    void buildInitialPopulation( population_init_type & pinit, unsigned int pN ) {
        if(! pinit.hasDistribution() ) return;
        
        size_t all_size = child_max_alleles( m_allele_space.size(), 0, pinit.getAlleleCount());

        // grow the child population accordingly
        m_child->grow( pN, all_size, m_trait_space.trait_count() );
        m_fit.resize( pN );

        // populate allele and trait space using configured generators
        generate_child_mutations( all_size );

        if( pinit.hasAlleleDistribution() ) {
            generateAlleleFrequency( pinit, pN, all_size );
        } else {
            generateGenotypeDistribution( pinit, pN, all_size );
        }
    }

    void generateGenotypeDistribution( population_init_type & pinit, unsigned int pN, unsigned int all_size ) {

        typedef typename population_init_type::genotype_dist_type::iterator iterator;

        std::vector< unsigned int > dist;
        dist.reserve( pN );
        unsigned all_idx = 0;
        for( iterator it = pinit.genotype_begin(); it != pinit.genotype_end(); it++ ) {
            clotho::genetics::genotype_distribution tmp = *it;
            assert( tmp.AA + tmp.AB + tmp.BA + tmp.BB == pN );

            dist.clear();

            for( unsigned int i = 0; i < tmp.AA; ++i ) {
                dist.push_back( 0 );
            }
            for( unsigned int i = 0; i < tmp.AB; ++i ) {
                dist.push_back( 1 );
            }
            for( unsigned int i = 0; i < tmp.BA; ++i ) {
                dist.push_back( 2 );
            }
            for( unsigned int i = 0; i < tmp.BB; ++i ) {
                dist.push_back( 3 );
            }

            std::shuffle( dist.begin(), dist.end(), *m_rand );

            unsigned int seq_idx = 0;
            for( std::vector< unsigned int >::iterator d_it = dist.begin(); d_it != dist.end(); d_it++ ) {
                if( *d_it == 0 || *d_it == 1 ) {
                    m_child->mutate( seq_idx, all_idx );
                }

                if( *d_it == 0 || *d_it == 3 ) {
                    m_child->mutate( seq_idx + 1, all_idx );
                }
                seq_idx += 2;
            }

            all_idx += 1;
        }
    }

    void generateAlleleFrequency( population_init_type & pinit, unsigned int pN, unsigned int all_size ) {
        double _pN = (double) m_child->haploid_genome_count();
        boost::random::uniform_int_distribution< unsigned int > seq_gen( 0, m_child->haploid_genome_count() - 1);

        unsigned int all_idx = 0;
        for( std::vector<double>::const_iterator it = pinit.allele_frequency_begin(); it != pinit.allele_frequency_end(); it++, ++all_idx ) {
            unsigned int N = (unsigned int) floor((*it) * _pN);

            for( unsigned int i = 0; i < N; ++i ) {
                unsigned int seq_idx = seq_gen( *m_rand );
                m_child->mutate(seq_idx, all_idx);
            }
        }

        bool allNeutral = m_trait_space.isAllNeutral();

        m_free_space.resetBuffers( m_child->getMaxBlocks() );
        free_buffer_type tbuf = m_free_space.getThreadBuffer(0);
        m_main_off_gen.reset( m_parent, m_child, tbuf, 0, pN, allNeutral );

        // generated population is post recombination and mutate
        // therefore only need to evaluate phenotype and fixed
        m_main_off_gen.phenotype();
        m_main_off_gen.fixed();
        
        m_fit( m_child );
        m_free_space.finalize(all_size);
    }

    void simulate( ) {
        std::swap( m_child, m_parent );     // use the current child population as the parent population for the next round

        // remove fixed alleles
        size_t free_count = removeFixedAlleles( m_parent );

        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = m_parent->getIndividualCount();
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
        m_fit.resize( pN );
        generate_child_mutations( pM );

        bool allNeutral = m_trait_space.isAllNeutral();

        thread_pool_type tpool( m_thread_count.m_tc );
        
        const size_t TC = m_thread_count.m_tc + 1;
        const size_t BATCH_SIZE = (pN / TC) + ((pN % TC > 0)? 1:0);

        m_free_space.resetBuffers( m_child->getMaxBlocks() );

        unsigned int off_idx = 0, j = 0;
        while( off_idx + BATCH_SIZE < pN ) {
            free_buffer_type tbuf = m_free_space.getThreadBuffer( j );
            m_worker_off_gens[ j ]->reset( m_parent, m_child, tbuf, off_idx, off_idx + BATCH_SIZE, allNeutral );
            off_idx += BATCH_SIZE;
            ++j;
        }

        tpool.post_list( m_worker_off_gens, m_thread_count.m_tc );

        if( off_idx < pN ) {
            assert( j < TC );
            free_buffer_type tbuf = m_free_space.getThreadBuffer( j );
            m_main_off_gen.reset( m_parent, m_child, tbuf, off_idx, pN, allNeutral );
            m_main_off_gen();
        }

        tpool.sync();

        timer_type fit_time;
        m_fit( m_child );
        fit_time.stop();

        clotho::utility::add_value_array( m_fitness_times, fit_time );

        m_free_space.finalize(all_size);

        recordAlleles();

        ++m_generation;
    }

    sequence_space_type * getChildPopulation() const {
        return m_child;
    }

    sequence_space_type * getParentPopulation() const {
        return m_parent;
    }

    void recordAlleles() {
        clotho::utility::add_value_array( free_sizes, m_free_space.free_size() );
        clotho::utility::add_value_array( var_sizes, m_free_space.variable_count() );
        clotho::utility::add_value_array( fixed_sizes, m_free_space.fixed_size() );
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

            boost::property_tree::ptree mut, xo, ph, fx;
            m_worker_off_gens[ i ]->record( xo, mut, ph, fx );
            
            log.put_child( "performance.mutate." + oss.str(), mut );
            log.put_child( "performance.crossover." + oss.str(), xo );
            log.put_child( "performance.fixed." + oss.str(), fx );
            log.put_child( "performance.phenotype." + oss.str(), ph );
        }

        boost::property_tree::ptree mut, xo, ph, fx;
        m_main_off_gen.record( xo, mut, ph, fx );
        
        log.put_child( "performance.mutate.main", mut );
        log.put_child( "performance.crossover.main", xo );
        log.put_child( "performance.fixed.main", fx );
        log.put_child( "performance.phenotype.main", ph );

        log.put_child( "performance.fitness", m_fitness_times);

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
        while( i < N ) {
            typename free_space_type::size_type all_idx = m_allele_space.size();

            if( it != end ) {
                all_idx = *it++;
            } else {
                m_allele_space.grow();
                m_trait_space.grow();
            }

            unsigned int seq_idx = seq_gen( *m_rand );

            assert( all_idx < m_child->getMaxAlleles() );

            m_mut_pool[ i ] = all_idx;
            m_mut_dist[ seq_idx + 1 ] += 1;

            allele_gen( m_allele_space, all_idx, m_generation );
            trait_gen(*m_rand, m_trait_space, all_idx );

            ++i;
        }

        // scan right to produce m_mut_pool relative index ranges
        for( unsigned int i = 2; i < m_mut_dist.size(); ++i ) {
            m_mut_dist[ i ] += m_mut_dist[ i - 1 ];
        }
    }

    void resetMutationEvents( unsigned int N ) {
        while( m_mut_pool.size() < N ) {
            m_mut_pool.push_back(0);
        }

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

    allele_type          m_allele_space, m_fixed;
    sequence_space_type  m_pop0, m_pop1;
    sequence_space_type  * m_parent, * m_child;

    trait_space_type        m_trait_space, m_fixed_traits;

//    selection_type          select_gen;
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

    offspring_generator_type m_main_off_gen;
    offspring_generator_type ** m_worker_off_gens;
//
    time_map m_fixed_times, m_mutate_times, m_xover_times, m_pheno_times;
    boost::property_tree::ptree free_sizes, var_sizes, fixed_sizes, m_fitness_times;
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
