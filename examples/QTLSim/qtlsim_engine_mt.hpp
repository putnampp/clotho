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

#include <boost/property_tree/ptree.hpp>

// including thread headers
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include "clotho/genetics/population_growth_toolkit.hpp"

#include "clotho/data_spaces/allele_space/allele_space.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"
#include "clotho/data_spaces/free_space/free_space.hpp"
#include "clotho/data_spaces/selection/selection.hpp"
#include "clotho/data_spaces/mutation/mutation_generator.hpp"
#include "clotho/data_spaces/mutation/mutation_allocator.hpp"

#include "clotho/data_spaces/crossover/crossover_task_list.hpp"

#include "clotho/data_spaces/phenotype_evaluator/evaluator.hpp"
#include "clotho/data_spaces/fitness/fitness.hpp"

#include "clotho/utility/state_object.hpp"

#include "thread_count_parameter.hpp"

template < class RNG, class RealType = double, class BlockType = unsigned long long >
class EngineMT {
public:
    typedef EngineMT< RNG >             self_type;

    typedef RealType                    position_type;
    typedef RealType                    weight_type;
    typedef weight_type *               phenotype_type;
    typedef BlockType                   block_type;

    typedef RNG                         random_engine_type;

    typedef clotho::genetics::qtl_allele_vectorized< position_type, weight_type >                   allele_type;

    typedef clotho::genetics::genetic_space< allele_type, block_type, __ALIGNMENT_TYPE__ > genetic_space_type;

    typedef clotho::genetics::mutation_allocator< random_engine_type, size_t >                      mutation_alloc_type;
    typedef clotho::genetics::MutationGenerator< random_engine_type, genetic_space_type >           mutation_type;
    typedef clotho::genetics::TraitWeightAccumulator< genetic_space_type >                          trait_accumulator_type;
    typedef clotho::genetics::FreeSpaceAnalyzer< genetic_space_type >                               free_space_type;
    typedef clotho::genetics::SelectionGenerator< random_engine_type, clotho::genetics::fitness_selection< genetic_space_type > >          selection_type;
//    typedef clotho::genetics::Crossover< random_engine_type, genetic_space_type >                   crossover_type;
    typedef clotho::genetics::linear_combination< trait_accumulator_type, phenotype_type >          trait_reduction_type;
    typedef clotho::genetics::phenotype_evaluator< trait_reduction_type >                           phenotype_eval_type;
    typedef clotho::genetics::Fitness< phenotype_eval_type >                                        fitness_type;

    typedef std::shared_ptr< ipopulation_growth_generator >                                         population_growth_generator_type;
    typedef std::shared_ptr< ipopulation_growth >                                                   population_growth_type;

    friend struct clotho::utility::state_getter< self_type >;

    EngineMT( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_parent( &m_pop0 )
        , m_child( &m_pop1 )
        , select_gen( rng, config )
        , mutate_gen( rng, config )
        , cross_gen( rng, config )
        , m_fit( config )
        , m_generation( 0 )
        , m_pop_growth()
        , m_mut_alloc( rng, config )
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

        m_pop0.getSequenceSpace().clear();
        m_pop1.getSequenceSpace().clear();

        size_t acount = m_pop0.getAlleleSpace().allele_count();
        m_pop0.getAlleleSpace().grow(acount, 1);
        m_pop1.getAlleleSpace().grow(acount, 1);

        thread_count_parameter tc_param( config );
        init(acount, tc_param.m_tc);
    }

    size_t getGeneration() const {
        return m_generation;
    }

    void init( size_t aN, int tc ) {
        size_t generation = m_generation++;
        size_t pN = m_parent->individual_count();
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, generation );
        }

        m_pop1.grow( pN, aN );

#ifdef USE_ROW_VECTOR
        m_pop1.getSequenceSpace().fill_empty();
        m_pop1.getSequenceSpace().finalize();
#endif // USE_ROW_VECTOR

        while( tc-- ) {
            threads.create_thread( boost::bind( &boost::asio::io_service::run, &m_service ) );
        }   
    }

    void simulate( ) {
        std::swap( m_child, m_parent );     // use the current child population as the parent population for the next round

        size_t generation = m_generation++;

        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = m_parent->individual_count();
        if( m_pop_growth ) {
            pN = m_pop_growth->operator()( pN, generation );
        }

        size_t pM = m_mut_alloc.allocate( 2 * pN );        // generate the number of new mutations

        BOOST_LOG_SEV( _log, boost::log::trivial::debug ) << "Generation " << generation << ": " << pN << " individuals; " << pM << " new alleles" << std::endl;

        select_gen.update( m_parent, pN );

        updateFixedAlleles( m_parent );                 // update the fixed alleles with those of parent population

        size_t all_size = child_max_alleles( m_parent->allele_count(), m_free_space.free_size(), pM );   // rescale allele space for child population given free space from parent population and new allele count (pM)
        BOOST_LOG_SEV( _log, boost::log::trivial::debug ) << "Rescaling child population to be: " << pN << " individuals x " << all_size << " alleles" << std::endl;

        m_child->grow( pN, all_size );                        // grow the child population accordingly

        m_child->inherit_alleles( m_parent, m_free_space.free_begin(), m_free_space.free_end() );

//        cross_gen.update( m_parent, select_gen.getMatePairs(), m_child );

        mutate_gen( m_child, pM );

        m_trait_accum.update( *m_child );
        m_pheno.update( m_child, m_trait_accum );
        m_fit.update( m_child, m_pheno );
    }

/*
    void dump( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree par, chi;

        m_parent->dump(par);
        l.put_child( "parent", par );

        m_child->dump(chi);
        l.put_child( "child", chi );

        boost::property_tree::ptree free;
        m_free_space.dump(free);
        l.put_child( "free", free );

        boost::property_tree::ptree tra;
        m_trait_accum.dump( tra );
        l.put_child( "trait", tra );
    
        boost::property_tree::ptree phe;
        m_pheno.dump(phe);
        l.put_child( "phenotype", phe);

        boost::property_tree::ptree fit;
        m_fit.update(fit);
        l.put_child( "fitness", fit );
    }*/

    genetic_space_type * getChildPopulation() const {
        return m_child;
    }

    genetic_space_type * getParentPopulation() const {
        return m_parent;
    }

    virtual ~EngineMT() {
        m_service.stop();

        m_threads.join_all();
    }

protected:

/**
 * estimate the maximum number of alleles in the child
 *
 * N_parent - number of alleles in the parent population
 * F_parent - number of free alleles in the parent population
 * M_child  - number of new alleles to be added the child population
 */
    size_t child_max_alleles( size_t N_parent, size_t F_parent, size_t M_child ) const {
//        std::cerr << "Parent alleles: " << N_parent << "; Free: " << F_parent << "; New Alleles: " << M_child << std::endl;

        if( F_parent >= M_child ) {
            // if there are more free alleles in the parent generation
            // than there are new alleles to be added to the child generation
            // then do not adjust scale of the allele space
            return N_parent;
        } else {
            return N_parent + (M_child - F_parent);
        }
    }

    void updateFixedAlleles( genetic_space_type * gs ) {
        m_free_space.update( *gs );               // analyze the parent population sequence space

        typedef typename free_space_type::iterator  fixed_iterator;

        fixed_iterator  fix_it = m_free_space.fixed_begin();
        fixed_iterator  fix_end = m_free_space.fixed_end();

        while( fix_it != fix_end ) {
            size_t  fixed_index = *fix_it++;

            gs->getSequenceSpace().flipColumn( fixed_index );

            m_fixed.push_back( gs->getAlleleSpace(), fixed_index );
        }

#ifdef DEBUGGING
        typedef typename free_space_type::iterator free_iterator;
        free_iterator fr_it = m_free_space.free_begin();
        free_iterator fr_end = m_free_space.free_end();

        while( fr_it != fr_end ) {
            size_t i = *fr_it++;

            if( !gs->getSequenceSpace().freeColumn( i ) ) {
                std::cerr << "Column is not actually free: " << i << " out of " << gs->getAlleleSpace().size() << std::endl;
                assert(false);
            }
        }
#endif // DEBUGGING
    }

    random_engine_type  * m_rand;
    genetic_space_type  m_pop0, m_pop1;
    genetic_space_type  * m_parent, * m_child;

    allele_type             m_fixed;

    trait_accumulator_type  m_trait_accum;
    phenotype_eval_type     m_pheno;
    free_space_type         m_free_space;

    selection_type          select_gen;
    mutation_type           mutate_gen;
    crossover_type          cross_gen;
    fitness_type            m_fit;

    size_t                  m_generation;

    population_growth_type  m_pop_growth;
    mutation_alloc_type     m_mut_alloc;

    // configuration parameter controlling the number of worker threads
    thread_count_parameter  m_multi;

    boost::asio::io_service m_service;
    boost::thread_grroup    m_threads;
};

namespace clotho {
namespace utility {

template < class RNG, class RealType, class BlockType >
struct state_getter< EngineMT< RNG, RealType, BlockType > > {
    typedef EngineMT< RNG, RealType, BlockType >           object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree tr;
        state_getter< typename object_type::trait_accumulator_type > tr_logger;
        tr_logger( tr, obj.m_trait_accum );

        boost::property_tree::ptree ph;
        state_getter< typename object_type::phenotype_eval_type > pheno_logger;
        pheno_logger( ph, obj.m_pheno );

        boost::property_tree::ptree fr;
        state_getter< typename object_type::free_space_type > free_logger;
        free_logger( fr, obj.m_free_space );

        boost::property_tree::ptree fx;
        state_getter< typename object_type::allele_type > all_logger;
        all_logger( fx, obj.m_fixed );

        boost::property_tree::ptree c_pop, p_pop;
        state_getter< typename object_type::genetic_space_type > pop_logger;
        pop_logger( p_pop, *(obj.m_parent) );
        pop_logger( c_pop, *(obj.m_child) );

        s.put_child( "traits", tr );
        s.put_child( "phenotypes", ph );
        s.put_child( "free_space", fr );
        s.put_child( "fixed_alleles", fx );

        //s.put_child( "parent", p_pop );
        s.put_child( "child", c_pop );
    }
};

}
}

#endif  // CLOTHO_SIM_ENGINE_MT_HPP_
