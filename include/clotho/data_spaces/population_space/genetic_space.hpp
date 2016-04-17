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
    typedef AlleleSpaceType                         allele_type;

    typedef BlockType                               block_type;
    typedef std::pair< size_t, size_t >             individual_genome_type;
    typedef std::vector< individual_genome_type >   population_genomes_type;

    typedef association_matrix< BlockType > association_type;
    typedef typename association_type::block_iterator   block_iterator;

    typedef typename association_type::row_iterator     sequence_iterator;

    genetic_space() {
        this->grow(1, 1);
    }

    allele_type & getAlleleSpace() {
        return m_alleles;
    }

    population_genomes_type & getIndividualSpace() {
        return m_genomes;
    }

    association_type & getSequenceSpace() {
        return m_assoc_data;
    }

    sequence_iterator    getSequenceAt( size_t idx ) const {
        return m_assoc_data.getRowAt( idx );
    }

    block_iterator  getBlockIterator() const {
        return m_assoc_data.raw_iterator();
    }

    size_t  allele_count() const {
        return m_alleles.size();
    }

    size_t  genome_count() const {
        return m_genomes.size();
    }

    size_t  sequence_count() const {
        return m_assoc_data.row_count();
    }

    individual_genome_type  getGenomeAt( size_t idx ) {
        assert( 0 <= idx && idx < m_genomes.size() );

        return m_genomes[ idx ];
    }

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

        if( genomes > m_genomes.size() ) {
            m_genomes.reserve( genomes );

            while( m_genomes.size() < genomes ) {
                m_genomes.push_back( std::make_pair( 0, 0) );
            }
        }

        alleles = m_alleles.grow( alleles );
    }

    allele_type           m_alleles;
    population_genomes_type     m_genomes;
    association_type            m_assoc_data;
};

}   // namespace clotho
}   // namespace genetics

#endif  // CLOTHO_GENETIC_SPACE_HPP_
