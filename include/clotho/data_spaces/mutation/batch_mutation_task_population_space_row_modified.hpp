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
#ifndef CLOTHO_BATCH_MUTATION_TASK_POPULATION_SPACE_ROW_MODIFIED_HPP_
#define CLOTHO_BATCH_MUTATION_TASK_POPULATION_SPACE_ROW_MODIFIED_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/population_space/population_space_row_modified.hpp"

#include "clotho/data_spaces/mutation/mutation_generators.hpp"
#include "clotho/data_spaces/allele_space/allele_generators.hpp"

#include "clotho/data_spaces/phenotype_evaluator/trait_space_generator.hpp"

#include <set>

namespace clotho {
namespace genetics {

template < class RNG, class BlockType, class WeightType, class AlleleSpaceType, class TraitSpaceType, class SizeType >
class BatchMutationTask< RNG, population_space_row_modified< BlockType, WeightType >, AlleleSpaceType, TraitSpaceType, SizeType > : public task {
public:

    typedef BatchMutationTask< RNG, population_space_row_modified< BlockType, WeightType >, AlleleSpaceType, TraitSpaceType, SizeType > self_type;

    typedef RNG random_engine_type;
    typedef population_space_row_modified< BlockType, WeightType > population_space_type;
    typedef AlleleSpaceType allele_space_type;
    typedef TraitSpaceType  trait_space_type;

    typedef SizeType    size_type;

    typedef std::vector< size_type > free_vector;
    typedef typename free_vector::iterator free_iterator;
    typedef typename free_vector::const_iterator const_free_iterator;

    typedef MutationGenerator2< random_engine_type, population_space_type >     mutation_generator_type;
    typedef typename mutation_generator_type::sequence_distribution_type        sequence_distribution_type;

    typedef AlleleGenerator< random_engine_type, allele_space_type >            allele_generator_type;

    typedef TraitSpaceGenerator2< trait_space_type >                            trait_generator_type;
    typedef typename trait_generator_type::parameter_type                       weight_parameter_type;

    BatchMutationTask( random_engine_type * rng, population_space_type * pop, allele_space_type * alleles, trait_space_type * traits, const_free_iterator first, const_free_iterator last, allele_generator_type & all_gen, unsigned int age, weight_parameter_type & wparam ) :
        m_rng( rng )
        , m_pop( pop )
        , m_alleles( alleles )
        , m_traits( traits )
        , m_free( first, last )
        , m_allele_gen( rng, all_gen )
        , m_trait_gen( wparam )
        , m_age( age )
    {}

    BatchMutationTask( const self_type & other ) :
        m_rng( other.m_rng )
        , m_pop( other.m_pop )
        , m_alleles( other.m_alleles )
        , m_traits( other.m_traits )
        , m_free( other.m_free )
        , m_allele_gen( other.m_allele_gen )
        , m_trait_gen( other.m_trait_gen )
        , m_age( other.m_age )
    {}

    void operator()() {
        sequence_distribution_type  seq_dist( 0, m_pop->haploid_genome_count() - 1 );
        const_free_iterator first = m_free.begin();

        std::set< int > seqs;
        while( first != m_free.end() ) {
            size_type all_idx = *first++;
            unsigned int seq_idx = seq_dist( *m_rng );

            m_pop->mutate( seq_idx, all_idx );

            m_allele_gen( *m_alleles, all_idx, m_age );
            m_trait_gen(  *m_rng, *m_traits, all_idx );

            seqs.insert( seq_idx );
        }

        m_pop->updateMutatedSequences( seqs.begin(), seqs.end() );
    }

    virtual ~BatchMutationTask() {}

protected:
    random_engine_type      * m_rng;
    population_space_type   * m_pop;
    allele_space_type       * m_alleles;
    trait_space_type        * m_traits;
    free_vector             m_free;

    allele_generator_type   m_allele_gen;
    trait_generator_type    m_trait_gen;

    unsigned int            m_age;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_BATCH_MUTATION_TASK_POPULATION_SPACE_ROW_MODIFIED_HPP_
