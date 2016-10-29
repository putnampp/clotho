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
#ifndef CLOTHO_BATCH_MUTATION_MT_GENERATOR_HPP_
#define CLOTHO_BATCH_MUTATION_MT_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/generators/weight_parameter.hpp"

#include "clotho/data_spaces/mutation/batch_mutation_tasks.hpp"
#include "clotho/data_spaces/mutation/mutation_allocator.hpp"


namespace clotho {
namespace genetics {

template < class RNG, class PopulationSpaceType, class AlleleSpaceType, class FreeSpaceType, class TraitSpaceType >
class BatchMutationMT {
public:
    typedef RNG random_engine_type;

    typedef PopulationSpaceType       space_type;
    typedef AlleleSpaceType allele_type;
    typedef FreeSpaceType   free_space_type;
    typedef TraitSpaceType  trait_space_type;

    typedef typename free_space_type::size_type size_type;

    typedef BatchMutationTask< RNG, PopulationSpaceType, AlleleSpaceType, TraitSpaceType, size_type >   task_type;

    typedef typename task_type::allele_generator_type                           allele_generator_type;
    typedef typename task_type::trait_generator_type::parameter_type            weight_parameter_type;

    typedef std::vector< size_type >                                            free_vector;

    typedef mutation_allocator< random_engine_type, size_type >                 allocator_type;

    BatchMutationMT( random_engine_type * rng, boost::property_tree::ptree & config) :
        m_rng( rng )
        , m_weight_param( config )
        , m_base_allele_gen( rng, config )
        , m_mut_allocator( rng, config )
    {}

    size_type generateNewMutation( unsigned int N_hap_genomes ) {
        return m_mut_allocator.allocate( N_hap_genomes );
    }

    template < class PoolType >
    void operator()( space_type * pop, allele_type * alleles, trait_space_type * traits, free_space_type * free_space, const unsigned int N, unsigned int age,  PoolType & pool ) {

        free_vector free_indices( free_space->free_begin(), free_space->free_end() );
        if( N > free_space->free_size() ) {
            const unsigned int new_alleles = (N - free_space->free_size());

            for( unsigned int b = 0; b < new_alleles; ++b ) {
                free_indices.push_back( alleles->size() );
                alleles->grow();
            }
        }

        batch_generate( pop, alleles, traits, free_indices, N, age, pool );
    }

    virtual ~BatchMutationMT() {}

protected:

    template < class PoolType >
    void batch_generate( space_type * pop, allele_type * alleles, trait_space_type * traits, const free_vector & free_indices, const unsigned int N, unsigned int age, PoolType & pool ) {
        const size_t TC = pool.pool_size() + 1;
        const size_t BATCH_SIZE = ( N / TC) + (( N % TC > 0) ? 1 : 0);

        unsigned int off_idx = 0, j = 0;
        while( off_idx + BATCH_SIZE < N ) {
            unsigned int off_end = off_idx + BATCH_SIZE;

            task_type x( pool.getRNG( j++ ), pop, alleles, traits, free_indices.begin() + off_idx, free_indices.begin() + off_end, m_base_allele_gen, age, m_weight_param );
            pool.post(x);

            off_idx = off_end;
        }

        if( off_idx < N ) {
            task_type t( m_rng, pop, alleles, traits, free_indices.begin() + off_idx, free_indices.begin() + N, m_base_allele_gen, age, m_weight_param );
            t();
        }

        pool.sync();
    }

    random_engine_type      * m_rng;

    weight_parameter_type   m_weight_param;

    allele_generator_type   m_base_allele_gen;

    allocator_type          m_mut_allocator;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_MUTATION_MT_GENERATOR_HPP_
