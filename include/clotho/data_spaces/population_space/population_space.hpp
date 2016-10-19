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
#ifndef CLOTHO_POPULATION_SPACE_HPP_
#define CLOTHO_POPULATION_SPACE_HPP_

#include <memory>
#include <vector>
#include <map>

#include "clotho/utility/bit_helper.hpp"

namespace clotho {
namespace genetics {

template < class WeightType >
class genetic_weight {
public:
    typedef genetic_weight< WeightType >            self_type;

    typedef WeightType                              weight_type;
    typedef std::vector< weight_type >              weight_vector;

    typedef typename weight_vector::iterator        weight_iterator;
    typedef typename weight_vector::const_iterator  const_weight_iterator;

    genetic_weight( unsigned int N = 1 ) {
        m_weights.reserve( N );
    }

    genetic_weight( const weight_vector & weights ) :
        m_weights( weights )
    {}

    genetic_weight( const self_type & other ) :
        m_weights( other.m_weights )
    {}

    const weight_vector & getWeights() const {
        return m_weights;
    }

    void update( const_weight_iterator first, const_weight_iterator last ) {
        if( m_weights.empty() ) {
            while( first != last ) {
                m_weights.push_back( *first++ );
            }
        } else {
            unsigned int i = 0;

            while( i < m_weights.size() && first != last ) {
                m_weights[ i++ ] += *first++;
            }
        }
    }

    weight_iterator begin_weight() {
        return m_weights.begin();
    }

    weight_iterator end_weight() {
        return m_weights.end();
    }

    const_weight_iterator begin_weight() const {
        return m_weights.begin();
    }

    const_weight_iterator end_weight() const {
        return m_weights.end();
    }

    virtual ~genetic_weight() {}

protected:
    weight_vector   m_weights;
};

template < class BlockType >
class genetic_sequence {
public:
    typedef genetic_sequence< BlockType >               self_type;

    typedef BlockType                                   block_type;
    typedef std::vector< block_type >                   sequence_vector;

    typedef typename sequence_vector::iterator          sequence_iterator;
    typedef typename sequence_vector::const_iterator    const_sequence_iterator;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    genetic_sequence( unsigned int N = 1 ) {
        m_data.reserve( N );
    }

    genetic_sequence( const sequence_vector & seq ) :
        m_data( seq )
    {}

    genetic_sequence( const self_type & other ) :
        m_data( other.m_data )
    {}

    static unsigned int calculateBlockCount( unsigned int M ) {
        return bit_helper_type::padded_block_count( M );
    }

    const sequence_vector & getSequence() const {
        return m_data;
    }

    void append_sequence( const block_type & b ) {
        m_data.push_back( b );
    }

    void set( unsigned int offset ) {
        unsigned int block_offset = rescale_to_offset( offset );

        m_data[ block_offset ] |= bit_helper_type::bit_offset( offset );
    }


    void unset( unsigned int offset ) {
        unsigned int block_offset = rescale_to_offset(offset);

        m_data[ block_offset ] &= ~( bit_helper_type::bit_offset( offset ) );
    }

    void remove_fixed( unsigned int bidx, block_type mask ) {
#ifdef DEBUGGING
        assert( m_data[ bidx ] & ~(mask) );       
#endif  // DEBUGGING
        if( bidx < m_data.size() )
            m_data[ bidx ] &= mask;
    }

    bool isSet( unsigned int idx ) {
        unsigned int b_idx = idx / bit_helper_type::BITS_PER_BLOCK;
        block_type mask = ((block_type) 1 << (idx % bit_helper_type::BITS_PER_BLOCK));
        return (m_data[ b_idx ] & mask );
    }

    bool isEmpty() const {
        return m_data.empty();
    }

    unsigned int sequence_block_length() const {
        return m_data.size();
    }

    sequence_iterator begin_sequence() {
        return m_data.begin();
    }

    sequence_iterator end_sequence() {
        return m_data.end();
    }

    const_sequence_iterator begin_sequence() const {
        return m_data.begin();
    }

    const_sequence_iterator end_sequence() const {
        return m_data.end();
    }

    virtual ~genetic_sequence() {}

protected:

    inline unsigned int rescale_to_offset( unsigned int offset ) {
        unsigned int block_offset = offset / bit_helper_type::BITS_PER_BLOCK;

        while( m_data.size() <= block_offset ) {
            //m_data.push_back( bit_helper_type::ALL_UNSET );
            m_data.push_back( 0 );
        }

        return block_offset;
    }

    sequence_vector m_data;
};

template < class BlockType, class WeightType >
class haploid_genome : public genetic_sequence< BlockType >, public genetic_weight< WeightType > {
public:
    typedef haploid_genome< BlockType, WeightType > self_type;

    typedef genetic_sequence< BlockType >           sequence_type;
    typedef genetic_weight< WeightType >            weight_type;

    haploid_genome( unsigned int N, unsigned int M) :
        sequence_type( N )
        , weight_type( M )
        , m_modifiable(true)
    {}

    haploid_genome( const self_type & other ) :
        sequence_type( other.getSequence() )
        , weight_type( other.getWeights() )
        , m_modifiable( true )
    {}

    bool isModifiable() const {
        return m_modifiable;
    }

    void finalize() {
        m_modifiable = false;
    }
   
    virtual ~haploid_genome() {}

protected:
    bool m_modifiable;
};

template < class BlockType, class WeightType >
class population_space {
public:

    typedef haploid_genome< BlockType, WeightType >         base_genome_type;

    typedef typename base_genome_type::sequence_type::sequence_vector     sequence_vector;

    typedef typename sequence_vector::iterator              sequence_iterator;
    typedef typename sequence_vector::const_iterator        const_sequence_iterator;

    typedef std::shared_ptr< base_genome_type >             genome_type;

    typedef std::pair< genome_type, genome_type >           individual_type;

    typedef std::vector< individual_type >                  population_type;

    typedef typename population_type::iterator              individual_iterator;
    typedef typename population_type::const_iterator        const_individual_iterator;

    typedef std::map< genome_type, unsigned int >           haploid_genomes;
//    typedef std::map< base_genome_type *, unsigned int >    haploid_genomes;
    typedef typename haploid_genomes::iterator              genome_iterator;
    typedef typename haploid_genomes::const_iterator        const_genome_iterator;

    typedef typename base_genome_type::weight_type::weight_type weight_type;

    population_space() :
        m_max_alleles(0)
        , m_max_blocks(0)
        , m_max_traits(0) 
    {}

    genome_type create_sequence( ) {
        return genome_type( new base_genome_type( m_max_blocks, m_max_traits ) );
    }

    typename base_genome_type::weight_type::weight_vector create_weight_vector( const weight_type default_value = 0.0  ) {
        return typename base_genome_type::weight_type::weight_vector( m_max_traits, default_value);
    }

    size_t individual_count() const {
        return m_pop.size();
    }

    size_t haploid_genome_count() const {
        return 2 * m_pop.size();
    }

    size_t haploid_genome_pointer_count() const {
        return m_genomes.size();
    }

    individual_type & getIndividual( size_t offset ) {
        return m_pop[ offset ];
    }

    void setIndividual( size_t offset, individual_type & ind ) {

        if( offset < m_pop.size() ) {
            m_pop[ offset ] = ind;
        } else {
            do {
                m_pop.push_back( ind );
            } while( m_pop.size() <= offset );
        }
    }

    void setIndividual( size_t offset, genome_type a, genome_type b ) {
        individual_type ind = std::make_pair( a, b );
        setIndividual( offset, ind );
    }

    individual_iterator begin_individual() {
        return m_pop.begin();
    }

    individual_iterator end_individual() {
        return m_pop.end();
    }

    const_individual_iterator begin_individual() const {
        return m_pop.begin();
    }

    const_individual_iterator end_individual() const {
        return m_pop.end();
    }

    genome_iterator begin_genomes() {
        return m_genomes.begin();
    }

    genome_iterator end_genomes() {
        return m_genomes.end();
    }

    const_genome_iterator begin_genomes() const {
        return m_genomes.begin();
    }

    const_genome_iterator end_genomes() const {
        return m_genomes.end();
    }

    void mutate( unsigned int seq_idx, unsigned int all_idx ) {
        unsigned int ind_idx = seq_idx / 2;

#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Mutate Individual: " << ind_idx << " [" << seq_idx << " vs " << m_pop.size() << "]";
#endif  // DEBUGGING
        assert( ind_idx < m_pop.size() );

        if( seq_idx % 2 ) {
            if( !m_pop[ ind_idx ].second ) {
                m_pop[ ind_idx ].second = create_sequence();
            } else if( !m_pop[ind_idx].second->isModifiable() ) {
                m_pop[ ind_idx ].second = genome_type( new base_genome_type(  *(m_pop[ ind_idx ].second) ) );
            }
            m_pop[ ind_idx ].second->set( all_idx );
        } else {
            if( !m_pop[ ind_idx ].first ) {
                m_pop[ ind_idx ].first = create_sequence();
            } else if( !m_pop[ind_idx].first->isModifiable() ) {
                m_pop[ ind_idx ].first = genome_type( new base_genome_type(  *(m_pop[ ind_idx ].first) ) );
            }
            m_pop[ ind_idx ].first->set( all_idx );
        }
    }

    void remove_fixed_allele( unsigned int all_idx ) {
        genome_iterator first = m_genomes.begin(), last = m_genomes.end();

        unsigned int bidx = all_idx / base_genome_type::sequence_type::bit_helper_type::BITS_PER_BLOCK;
        typename base_genome_type::sequence_type::block_type mask = ~(base_genome_type::sequence_type::bit_helper_type::bit_offset( all_idx ));
        while( first != last ) {
            if(first->first) first->first->remove_fixed( bidx, mask );
            ++first;
        }
    }

    bool freeColumn( unsigned int idx ) {
        genome_iterator first = m_genomes.begin(), last = m_genomes.end();

        bool res = true;
        while( res && first != last ) {
            if( !first->first ) { 
                res = first->first->isSet( idx );
            }
            ++first;
        }

        return res;
    }

    void grow( unsigned int I, unsigned int A ) {
        clear();

        std::cerr << "Reserving: " << I << std::endl;
        m_pop.reserve( I );

        // fill population with empty individuals
        while( m_pop.size() < I ) {
//            m_pop.push_back( std::make_pair( genome_type(), genome_type() ) );
            m_pop.push_back( std::make_pair( create_sequence(), create_sequence()) );
        }

        setMaxAlleles( A );
    }
/*
    void finalize() {
        individual_iterator it = m_pop.begin();

        while( it != m_pop.end() ) {
            if( it->first ) it->first->finalize();
            if( it->second ) it->second->finalize();

            genome_iterator res = m_genomes.find( it->first.get() );
            if( res == m_genomes.end() ) {
                m_genomes.insert( std::make_pair( it->first.get(), 1 ) );
            } else {
                res->second += 1;
            }

            res = m_genomes.find( it->second.get() );
            if( res == m_genomes.end() ) {
                m_genomes.insert( std::make_pair( it->second.get(), 1 ) );
            } else {
                res->second += 1;
            }

            ++it;
        }
    }
*/
    void finalize() {
        individual_iterator it = m_pop.begin();

        while( it != m_pop.end() ) {
            if( it->first ) it->first->finalize();
            if( it->second ) it->second->finalize();

            genome_iterator res = m_genomes.find( it->first );
            if( res == m_genomes.end() ) {
                m_genomes.insert( std::make_pair( it->first, 1 ) );
            } else {
                res->second += 1;
            }

            res = m_genomes.find( it->second );
            if( res == m_genomes.end() ) {
                m_genomes.insert( std::make_pair( it->second, 1 ) );
            } else {
                res->second += 1;
            }

            ++it;
        }
    }

    void setMaxAlleles( unsigned int M ) {
        m_max_alleles = M;
        m_max_blocks = base_genome_type::calculateBlockCount( M );
    }

    void setMaxTraits( unsigned int T ) {
        m_max_traits = T;
    }

    unsigned int getMaxAlleles( ) const {
        return m_max_alleles;
    }

    unsigned int getMaxBlocks( ) const {
        return m_max_blocks;
    }

    unsigned int getMaxTraits() const {
        return m_max_traits;
    }

    void clear() {
        m_pop.clear();
        m_genomes.clear();
    }

    virtual ~population_space() {}

protected:
    population_type     m_pop;
    haploid_genomes     m_genomes;

    unsigned int        m_max_alleles, m_max_blocks, m_max_traits;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POPULATION_SPACE_HPP_
