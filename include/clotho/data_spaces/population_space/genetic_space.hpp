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
#ifndef CLOTHO_GENETIC_SPACE_HPP_
#define CLOTHO_GENETIC_SPACE_HPP_

#include "clotho/data_spaces/growable2D.hpp"
#include "clotho/data_spaces/association_matrix.hpp"
#include <vector>

#include "clotho/utility/state_object.hpp"

namespace clotho {
namespace genetics {

template < class AlleleSpaceType, class BlockType, class OrganizationType >
class genetic_space : public growable2D {
public:
    typedef genetic_space< AlleleSpaceType, BlockType, OrganizationType >     self_type;
    typedef AlleleSpaceType                                 allele_type;

    typedef BlockType                                       block_type;
    typedef size_t                                          sequence_id_type;
    typedef size_t                                          individual_id_type;
    
    typedef std::pair< sequence_id_type, sequence_id_type > individual_genome_type;
    typedef std::vector< individual_genome_type >           population_type;

    typedef typename population_type::iterator          individual_iterator;
    typedef typename population_type::const_iterator    const_individual_iterator;

    typedef association_matrix< BlockType, OrganizationType >   association_type;
//    typedef typename association_type::block_iterator           block_iterator;

    typedef typename association_type::raw_block_pointer         sequence_iterator;
    typedef typename association_type::raw_block_pointer    genome_iterator;

    typedef double                                          fitness_score_type;
    typedef std::vector< fitness_score_type >               fitness_scores;
    typedef typename fitness_scores::iterator               fitness_iterator;
    typedef typename fitness_scores::const_iterator         const_fitness_iterator;

    genetic_space() {
        this->grow(1, 1);
    }

    allele_type & getAlleleSpace() {
        return m_alleles;
    }

    association_type & getSequenceSpace() {
        return m_assoc_data;
    }

    fitness_scores & getFitnessSpace() {
        return m_fitness;
    }

//    sequence_iterator    getSequenceAt( size_t idx ) const {
//        return m_assoc_data.getRowAt( idx );
//    }
//
//    block_iterator  getBlockIterator() const {
//        return m_assoc_data.raw_iterator();
//    }

    template < class FreeIterator >
    void inherit_alleles( self_type * parent_pop, FreeIterator first, FreeIterator last ) {

        m_alleles.inherit( parent_pop->m_alleles );
        m_alleles.updateFreeSpace( first, last );
    }

    sequence_iterator begin_sequence( size_t idx ) const {
        return m_assoc_data.begin_row( idx );
    }

    sequence_iterator end_sequence( size_t idx ) const {
        return m_assoc_data.end_row( idx );
    }

    genome_iterator begin_genome( individual_id_type idx ) const {
        return m_assoc_data.begin_row(2 * idx);
    }

    genome_iterator end_genome( individual_id_type idx ) const {
        return m_assoc_data.end_row(2 * idx);
    }

    size_t  allele_count() const {
        return m_alleles.size();
    }

    size_t  individual_count() const {
        return m_assoc_data.row_count() / 2;
    }

    size_t  genome_count() const {
        return m_assoc_data.row_count() / 2;
    }

    size_t  sequence_count() const {
        return m_assoc_data.row_count();
    }

    individual_genome_type  getIndividualAt( individual_id_type idx ) {
        assert( 0 <= idx && idx < individual_count() );

        return std::make_pair( 2 * idx, 2 * idx + 1);
    }

    void setFitnessAt( size_t idx, fitness_score_type f ) {
        assert( 0 <= idx && idx < individual_count() );

        m_fitness[ idx ] = f;
    }

    fitness_iterator    fitness_begin() {
        return m_fitness.begin();
    }

    fitness_iterator    fitness_end() {
        return m_fitness.end();
    }

    const_fitness_iterator    fitness_begin() const {
        return m_fitness.begin();
    }

    const_fitness_iterator    fitness_end() const {
        return m_fitness.end();
    }

    virtual size_t  grow( size_t individuals, size_t alleles) {
        this->resize( individuals, alleles );
        return 0;
    }

    virtual ~genetic_space() {}
protected:

    void resize( size_t individuals, size_t alleles ) {
        // 1 association row represents 1 'chromosome'
        // 2 chromosomes per genome
        //
        size_t seqs = 2 * individuals;
        m_assoc_data.grow( seqs, alleles );

        m_alleles.grow( alleles );

        m_fitness.resize( individuals );
    }

    allele_type                 m_alleles;
    association_type            m_assoc_data;

    fitness_scores              m_fitness;
};

}   // namespace clotho
}   // namespace genetics

namespace clotho {
namespace utility {


template < class AlleleSpaceType, class BlockType, class OrganizationType >
struct state_getter< clotho::genetics::genetic_space< AlleleSpaceType, BlockType, OrganizationType > > {
    typedef clotho::genetics::genetic_space< AlleleSpaceType, BlockType, OrganizationType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {

        boost::property_tree::ptree al;
        state_getter< typename object_type::allele_type > a_log;
        a_log( al, obj.getAlleleSpace() );

        boost::property_tree::ptree ft;
        clotho::utility::add_value_array( ft, obj.fitness_begin(), obj.fitness_end() );

        s.put_child( "alleles" , al );
        s.put_child( "fitness", ft );
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_GENETIC_SPACE_HPP_
