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

namespace clotho {
namespace genetics {

template < class AlleleSpaceType, class BlockType = unsigned long long >
class genetic_space : public growable2D {
public:
    typedef AlleleSpaceType                                 allele_type;

    typedef BlockType                                       block_type;
    typedef size_t                                          sequence_id_type;
    typedef size_t                                          individual_id_type;
    
    typedef std::pair< sequence_id_type, sequence_id_type > individual_genome_type;
    typedef std::vector< individual_genome_type >           population_type;

    typedef typename population_type::iterator          individual_iterator;
    typedef typename population_type::const_iterator    const_individual_iterator;

    typedef association_matrix< BlockType >                 association_type;
    typedef typename association_type::block_iterator       block_iterator;

    typedef typename association_type::row_iterator         sequence_iterator;
    typedef typename association_type::row_pair_iterator    genome_iterator;

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

//    population_type & getIndividualSpace() {
//        return m_population;
//    }

    association_type & getSequenceSpace() {
        return m_assoc_data;
    }

    fitness_scores & getFitnessSpace() {
        return m_fitness;
    }

    sequence_iterator    getSequenceAt( size_t idx ) const {
        return m_assoc_data.getRowAt( idx );
    }

    block_iterator  getBlockIterator() const {
        return m_assoc_data.raw_iterator();
    }

    genome_iterator getGenomeAt( individual_id_type idx ) const {
        return m_assoc_data.getRowPairAt(2 * idx );
    }

    size_t  allele_count() const {
        return m_alleles.size();
    }

    size_t  individual_count() const {
        return m_assoc_data.row_count() / 2;
//        return m_population.size();
    }

    size_t  genome_count() const {
//        return m_population.size();
        return m_assoc_data.row_count() / 2;
    }

    size_t  sequence_count() const {
        return m_assoc_data.row_count();
    }

    individual_genome_type  getIndividualAt( individual_id_type idx ) {
        assert( 0 <= idx && idx < individual_count() );

        return std::make_pair( 2 * idx, 2 * idx + 1);
//        return m_population[ idx ];
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

//    individual_iterator         individual_begin() {
//        return m_population.begin();
//    }
//
//    individual_iterator         individual_end() {
//        return m_population.end();
//    }
//
//    const_individual_iterator         individual_begin() const {
//        return m_population.begin();
//    }
//
//    const_individual_iterator         individual_end() const {
//        return m_population.end();
//    }
//
    virtual size_t  grow( size_t genomes, size_t alleles) {
        this->resize( genomes, alleles );
        return 0;
    }

    virtual ~genetic_space() {}
protected:

    void resize( size_t genomes, size_t alleles ) {
        // 1 association row represents 1 'chromosome'
        // 2 chromosomes per genome
        //
        m_assoc_data.grow( 2 * genomes, alleles );

        m_alleles.grow( alleles );

//        if( genomes > m_population.size() ) {
//            m_population.reserve( genomes );
//
//            while( m_population.size() < genomes ) {
//                m_population.push_back( std::make_pair( 0, 0) );
//            }
//        }

        alleles = m_alleles.grow( alleles );
    }

    allele_type                 m_alleles;
//    population_type     m_population;
    association_type            m_assoc_data;

    fitness_scores              m_fitness;
};

}   // namespace clotho
}   // namespace genetics

#endif  // CLOTHO_GENETIC_SPACE_HPP_
