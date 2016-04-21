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
#ifndef CLOTHO_SIM_ENGINE_HPP_
#define CLOTHO_SIM_ENGINE_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/data_spaces/allele_space/qtl_allele.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"
#include "clotho/data_spaces/free_space/free_space.hpp"
#include "clotho/data_spaces/selection/selection_generator.hpp"
#include "clotho/data_spaces/crossover/crossover.hpp"

template < class RNG >
class Engine {
public:
    typedef double                      position_type;
    typedef double                      weight_type;
    typedef unsigned long long          block_type;

    typedef RNG                                                             random_engine_type;

    typedef clotho::genetics::qtl_allele< position_type, weight_type >      allele_type;
    typedef clotho::genetics::genetic_space< allele_type, block_type >      genetic_space_type;

    typedef clotho::genetics::TraitWeightAccumulator< genetic_space_type >  trait_accumulator_type;
    typedef clotho::genetics::FreeSpaceAnalyzer< genetic_space_type >       free_space_type;
    typedef clotho::genetics::SelectionGenerator< genetic_space_type >      selection_type;
    typedef clotho::genetics::Crossover< random_engine_type, genetic_space_type >               crossover_type;

    Engine( random_engine_type * rng, boost::property_tree::ptree & config ) : 
        m_parent( &m_pop0 )
        , m_child( &m_pop1 )
        , select_gen( config )
        , cross_gen( rng, config )
    {}

    void simulate() {
        std::swap( m_child, m_parent );     // use the current child population as the parent population for the next round
        // at the start of each simulate round, m_fit has already been updated from the previous
        // round with the fitness of the "then child/now parent" popualtions fitness
        //
        size_t pN = select_gen.update( m_parent, m_fit );    // generate the number of children
        size_t pM = mutate_gen.generateSize( m_parent );     // generate the number of new mutations

        updateFixedAlleles( m_parent );                 // update the fixed alleles with those of parent population

        pM = child_max_alleles( m_parent->allele_count(), m_free_space.free_size(), pM );

        m_child->grow( pN, pM );                        // grow the child population accordingly

        m_child->updateFreeSpace( m_free.free_begin(), m_free.free_end() );

        cross_gen.update( m_parent, select_gen.getMatePairs(), m_child );

        mut_gen.update( m_child, m_parent );

        m_trait_accum.update( *m_child );
        m_pheno.update( *m_child, m_trait_accum );
        m_fit.update( m_pheno );
    }

    genetic_space_type * getChildPopulation() const {
        return m_child;
    }

    genetic_space_type * getParentPopulation() const {
        return m_parent;
    }

    virtual ~Engine() {}

protected:

/**
 * estimate the maximum number of alleles in the child
 *
 * N_parent - number of alleles in the parent population
 * F_parent - number of free alleles in the parent population
 * M_child  - number of new alleles to be added the child population
 */
    size_t child_max_alleles( size_t N_parent, size_t F_parent, size_t M_child ) const {
        return (N_parent - F_parent) + M_child;
    }

    void updateFixedAlleles( genetic_space_type * gs ) {
        m_free_space.update( *gs );               // analyze the parent population
    }

    genetic_space_type m_pop0, m_pop1;
    genetic_space_type * m_parent, * m_child;

    allele_type             m_fixed;

    trait_accumulator_type  m_trait_accum;
    free_space_type         m_free_space;

    selection_type          select_gen;
    mutation_type           mutate_gen;
    crossover_type          cross_gen;
};

#endif  // CLOTHO_SIM_ENGINE_HPP_
