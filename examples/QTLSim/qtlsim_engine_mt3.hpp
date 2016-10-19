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
#ifndef CLOTHO_SIM_ENGINE_MT_HPP_
#define CLOTHO_SIM_ENGINE_MT_HPP_

#ifdef DEBUG_MODE
#define DEBUGGING 0
#endif  // DEBUG_MODE

#include "qtlsim_logger.hpp"

#include <boost/property_tree/ptree.hpp>

#include "clotho/genetics/population_growth_toolkit.hpp"

#include "clotho/data_spaces/allele_space/allele_space_vector.hpp"
#include "clotho/data_spaces/allele_space/allele_generator_vector.hpp"

#include "clotho/data_spaces/phenotype_evaluator/trait_space_vector.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_space_generator.hpp"

#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"
#include "clotho/data_spaces/free_space/free_space_mts.hpp"

#ifdef USE_BATCH_JOBS
#include "clotho/data_spaces/crossover/batch_crossover_mt.hpp"
#define CROSSOVER_TYPE clotho::genetics::BatchCrossoverMT

#include "clotho/data_spaces/phenotype_evaluator/batch_phenotype_mts.hpp"
#define PHENOTYPE_TYPE clotho::genetics::BatchPhenotypeMT

#else
#include "clotho/data_spaces/crossover/crossover_mt.hpp"
#define CROSSOVER_TYPE clotho::genetics::CrossoverMT

#include "clotho/data_spaces/phenotype_evaluator/phenotype_mt.hpp"
#define PHENOTYPE_TYPE clotho::genetics::PhenotypeMT
#endif  // USE_BATCH_JOBS

#include "clotho/data_spaces/population_space/population_spaces.hpp"
#include "clotho/data_spaces/selection/selection.hpp"

#include "clotho/data_spaces/mutation/mutation_generators.hpp"
#include "clotho/data_spaces/mutation/mutation_allocator.hpp"

#include "clotho/data_spaces/fitness/general_fitness.hpp"

#include "clotho/utility/state_object.hpp"

#include "clotho/data_spaces/task/thread_pool.hpp"

template < class RNG, class RealType = double, class BlockType = unsigned long long, class SizeType = unsigned int >
class EngineMT2 {
public:
    typedef EngineMT2< RNG >             self_type;

    typedef RealType                    position_type;
    typedef RealType                    weight_type;
    typedef weight_type *               phenotype_type;
    typedef BlockType                   block_type;

    typedef RNG                         random_engine_type;

    typedef SizeType                    size_type;

    typedef clotho::genetics::thread_pool< RNG >                                                    thread_pool_type;
    typedef clotho::genetics::AlleleSpace< position_type, size_type >                               allele_type;

//    typedef clotho::genetics::population_space_columnar< block_type, weight_type >                  sequence_space_type;
    typedef clotho::genetics::population_space_row< block_type, weight_type >                  sequence_space_type;
    typedef clotho::genetics::trait_space_vector< weight_type >                                     trait_space_type;
    typedef clotho::genetics::FreeSpaceAnalyzerMT< sequence_space_type, size_type >                 free_space_type;

    typedef clotho::genetics::mutation_allocator< random_engine_type, size_type >                 mutation_alloc_type;

    typedef clotho::genetics::MutationGenerator2< random_engine_type, sequence_space_type >       mutation_type;

    typedef clotho::genetics::AlleleGenerator< random_engine_type, allele_type >                  allele_generator_type;
    typedef clotho::genetics::TraitSpaceGenerator< random_engine_type, trait_space_type >         trait_generator_type;

    typedef clotho::genetics::GeneralFitness                                                      fitness_type;
    typedef clotho::genetics::SelectionGenerator< random_engine_type, clotho::genetics::fitness_selection< fitness_type > >          selection_type;
    typedef CROSSOVER_TYPE< random_engine_type, sequence_space_type, allele_type >                crossover_type;

    typedef PHENOTYPE_TYPE< sequence_space_type, trait_space_type >                               phenotype_eval_type;

    typedef std::shared_ptr< ipopulation_growth_generator >                                         population_growth_generator_type;
    typedef std::shared_ptr< ipopulation_growth >                                                   population_growth_type;

    friend struct clotho::utility::state_getter< self_type >;

    EngineMT2( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_parent( &m_pop0 )
        , m_child( &m_pop1 )
        , m_trait_space( config )
        , m_fixed_traits( config )
        , m_thread_pool( rng, config )
        , m_free_space( )
        , select_gen( rng, config )
        , mutate_gen( rng, config )
        , cross_gen( rng, config )
        , m_fit( config )
        , m_generation( 0 )
        , m_pop_growth()
        , m_mut_alloc( rng, config )
        , trait_gen( rng, config )
        , allele_gen( rng, config )
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

//        m_pop0.getSequenceSpace().clear();
//        m_pop1.getSequenceSpace().clear();

//        size_t acount = m_pop0.getAlleleSpace().allele_count();
//        m_pop0.getAlleleSpace().grow(acount, m_pheno.trait_count());
//        m_pop1.getAlleleSpace().grow(acount, m_pheno.trait_count());

        init(0);
    }

    size_t getGeneration() const {
        return m_generation;
    }

    void init( size_t aN ) {
//        size_t pN = m_parent->individual_count();
        size_t pN = 0;
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        m_pop1.grow( pN, aN, m_trait_space.trait_count() );

#ifdef USE_ROW_VECTOR
        m_pop1.getSequenceSpace().fill_empty();
        m_pop1.getSequenceSpace().finalize();
#endif // USE_ROW_VECTOR
        ++m_generation;
    }

    void simulate( ) {
        std::swap( m_child, m_parent );     // use the current child population as the parent population for the next round

        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = select_gen.individual_count();
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, m_generation );
        }

        size_type pM = m_mut_alloc.allocate( 2 * pN );        // generate the number of new mutations
        size_type free_count = updateFixedAlleles( m_parent );    // update the fixed alleles with those of parent population
        size_t all_size = child_max_alleles( m_allele_space.size(), free_count, pM );   // rescale allele space for child population given free space from parent population and new allele count (pM)

#ifdef DEBUGGING 
        BOOST_LOG_TRIVIAL( debug ) << "Generation " << m_generation << ": " << pN << " individuals; " << pM << " new alleles";
        BOOST_LOG_TRIVIAL( debug ) << "Free space: " << free_count << "; alleles: " << m_allele_space.size();
        BOOST_LOG_TRIVIAL( debug ) << "Rescaling child population to be: " << pN << " individuals x " << all_size << " alleles";
        std::cerr << "Generation " << m_generation << ": " << pN << " individuals; " << pM << " new alleles" << std::endl;
        std::cerr << "Rescaling child population to be: " << pN << " individuals x " << all_size << " alleles" << std::endl;
#endif // DEBUGGING

        m_child->grow( pN, all_size, m_trait_space.trait_count() );               // grow the child population accordingly

        select_gen.update( m_fit, pN );

        cross_gen( select_gen, m_parent, m_child, &m_allele_space, m_thread_pool );

        generate_child_mutations( pM );

        if( !m_allele_space.isAllNeutral() ) {
            m_pheno( m_child, &m_trait_space, m_thread_pool );
        } else {
            m_pheno.constant_phenotype( m_child, &m_trait_space );
        }
        m_fit( m_pheno );

        ++m_generation;
    }

    sequence_space_type * getChildPopulation() const {
        return m_child;
    }

    sequence_space_type * getParentPopulation() const {
        return m_parent;
    }

    allele_type * getAlleleSpace() {
        return &m_allele_space;
    }

    virtual ~EngineMT2() { }

protected:

    void generate_child_mutations( unsigned int N ) {
//        std::cerr << "Child population size: " << m_child->haploid_genome_count() << std::endl;
        typename mutation_type::sequence_distribution_type seq_gen( 0, m_child->haploid_genome_count() - 1);

        typename free_space_type::base_type::iterator it = m_free_space.free_begin(), end = m_free_space.free_end();
        while( N && it != end ) {
            typename free_space_type::size_type all_idx = *it++;
            unsigned int seq_idx = seq_gen( *m_rand );

            mutate_gen( m_child, seq_idx, all_idx );
            allele_gen( m_allele_space, all_idx, m_generation );
            trait_gen( m_trait_space, all_idx );
            --N;
        }

        while( N ) {
            typename free_space_type::size_type all_idx = m_allele_space.size();
            unsigned int seq_idx = seq_gen( *m_rand );

            assert( all_idx < m_child->getMaxAlleles() );

            mutate_gen( m_child, seq_idx, all_idx );
            allele_gen( m_allele_space, all_idx, m_generation );
            trait_gen( m_trait_space, all_idx );
            --N;
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

    size_type updateFixedAlleles( sequence_space_type * ss ) {
        m_free_space( ss, m_thread_pool );               // analyze the parent population sequence space

        typedef typename free_space_type::iterator  fixed_iterator;
        typedef typename trait_space_type::iterator trait_iterator;


//        std::cerr << "Fixed count: " << m_free_space.fixed_size() << std::endl;

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

    thread_pool_type        m_thread_pool;
    phenotype_eval_type     m_pheno;
    free_space_type         m_free_space;

    selection_type          select_gen;
    mutation_type           mutate_gen;
    crossover_type          cross_gen;
    fitness_type            m_fit;

    size_t                  m_generation;

    population_growth_type  m_pop_growth;
    mutation_alloc_type     m_mut_alloc;

    trait_generator_type    trait_gen;
    allele_generator_type   allele_gen;
};

namespace clotho {
namespace utility {

template < class RNG, class RealType, class BlockType, class SizeType >
struct state_getter< EngineMT2< RNG, RealType, BlockType, SizeType > > {
    typedef EngineMT2< RNG, RealType, BlockType, SizeType >           object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree tr;
        state_getter< typename object_type::trait_space_type > tr_logger;
        tr_logger( tr, obj.m_trait_space );

        boost::property_tree::ptree ph;
        state_getter< typename object_type::phenotype_eval_type > pheno_logger;
        pheno_logger( ph, obj.m_pheno );

        boost::property_tree::ptree fr;
        state_getter< typename object_type::free_space_type > free_logger;
        free_logger( fr, obj.m_free_space );

        boost::property_tree::ptree fx, alls;
        state_getter< typename object_type::allele_type > all_logger;
        all_logger( fx, obj.m_fixed );
        all_logger( alls, obj.m_allele_space );

//        boost::property_tree::ptree c_pop;
//        state_getter< typename object_type::sequence_space_type > pop_logger;
//        pop_logger( c_pop, *(obj.m_child) );

        s.put_child( "phenotypes", ph );
        s.put_child( "free_space", fr );
        s.put_child( "allele_space", alls );
        s.put_child( "trait_space", tr );
        s.put_child( "fixed_alleles", fx );

//        s.put_child( "child", c_pop );
    }
};

}
}

#endif  // CLOTHO_SIM_ENGINE_MT_HPP_
