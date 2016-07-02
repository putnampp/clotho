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
#ifndef CLOTHO_ROW_GROUPED_CROSSOVER_METHOD_HPP_
#define CLOTHO_ROW_GROUPED_CROSSOVER_METHOD_HPP_

#include "clotho/data_spaces/crossover/crossover_method_def.hpp"
#include "clotho/data_spaces/association_matrix/row_grouped_association_matrix.hpp"

#include "clotho/data_spaces/generators/small_crossover_event_generator.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class AlleleSpaceType, class BlockType >
class crossover_method< RNG, AlleleSpaceType, association_matrix< BlockType, row_grouped< 2 > > > {
public:
    typedef RNG                                                 random_engine_type;
    typedef AlleleSpaceType                                     allele_type;
    typedef association_matrix< BlockType, row_grouped< 2 > >     sequence_space_type;

    typedef typename sequence_space_type::raw_block_pointer     sequence_iterator;
    typedef typename sequence_space_type::raw_block_pointer     genome_iterator;

    typedef typename allele_type::position_type                 position_type;

    typedef crossover_event_generator< RNG, position_type >     event_generator_type;

    typedef typename sequence_space_type::block_type            block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    crossover_method( RNG * rng, boost::property_tree::ptree & config ) :
        m_event_gen( rng, config )
    {}

    void initAlleles( allele_type & alls ) {
        m_event_gen.update( alls.getPositions(), alls.size());
    }

    void crossover( genome_iterator start, genome_iterator end, sequence_iterator s, size_t p_step, size_t s_step ) {

        size_t N = m_event_gen.generate();

        if( N != 0 ) {
            size_t j = 0;
            while( start != end ) {
                block_type first = *start++;
                block_type second = *start++;

                block_type hets = first ^ second;
                block_type sec  = bit_helper_type::ALL_UNSET;   // mask state from second strand

                while( hets ) {
                    size_t idx = bit_walker_type::unset_next_index( hets ) + j;

                    if( m_event_gen( idx ) ) {
                        sec |= bit_helper_type::bit_offset( idx );
                    }
                }

                *s = ((first & ~sec) | (second & sec));

                j += bit_helper_type::BITS_PER_BLOCK;

                s += row_grouped< 2 >::GROUP_SIZE;
            }
        } else {
            if( m_event_gen.getBaseSequence() != 0 ) {

                ++start;
                ++end;
            }

            size_t j = 0;
            while( start != end ) {
                *s = *start;

                start += row_grouped< 2 >::GROUP_SIZE;
                s += row_grouped< 2 >::GROUP_SIZE;

                j += bit_helper_type::BITS_PER_BLOCK;
            }
        }
    }

    virtual ~crossover_method() {}

protected:
    event_generator_type m_event_gen;
};

template < class RNG, class AlleleSpaceType, class BlockType >
class crossover_method< RNG, AlleleSpaceType, association_matrix< BlockType, row_grouped< 1 > > > {
public:
    typedef RNG                                                 random_engine_type;
    typedef AlleleSpaceType                                     allele_type;
    typedef association_matrix< BlockType, row_grouped< 1 > >   sequence_space_type;

    typedef typename sequence_space_type::row_vector            row_vector;
    typedef typename row_vector::raw_pointer                    raw_pointer;

    typedef typename allele_type::position_type                 position_type;

    typedef small_crossover_event_generator< RNG, position_type >     event_generator_type;

    typedef typename sequence_space_type::block_type            block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    crossover_method( RNG * rng, boost::property_tree::ptree & config ) :
        m_event_gen( rng, config )
        , m_size(0)
        , m_count(0)
        , m_buffer( NULL )
    {
#ifdef DEBUGGING
        std::cerr << "Using specialized crossover method" << std::endl;
#endif  // DEBUGGING
    }

    void initAlleles( allele_type & alls ) {
        m_event_gen.update( alls.getPositions(), alls.size() );
        m_size = bit_helper_type::padded_block_count( alls.size() );
        m_count = alls.size();

        if( m_buffer != NULL ) {
            delete [] m_buffer;
            m_buffer = NULL;
        }
    }

    row_vector crossover( row_vector & a, row_vector & b ) {

        size_t N = m_event_gen.generate();

        if( N != 0 ) {
            if( a == b ) {
                return row_vector(a);
            }

            if( m_buffer == NULL ) {
                m_buffer = new block_type[ m_size ];
            }

            size_t N = a.m_size, M = b.m_size;

#ifdef DEBUGGING
            assert(N <= m_size && M <= m_size);
            assert( a.m_readonly );
            assert( b.m_readonly );
#endif  // DEBUGGING

            raw_pointer a_ptr = a.get(), b_ptr = b.get();
            raw_pointer c_ptr = m_buffer;

            bool match_first = true, match_second = true;
            size_t i = 0, j = 0;
            while( i < m_size ) {
                block_type first = (( i < N ) ? a_ptr[ i ] : bit_helper_type::ALL_UNSET);
                block_type second = (( i < M ) ? b_ptr[ i ] : bit_helper_type::ALL_UNSET);

                block_type hets = first ^ second;
                block_type sec  = bit_helper_type::ALL_UNSET;   // mask state from second strand

                while( hets ) {
                    size_t idx = bit_walker_type::unset_next_index( hets ) + j;

#ifdef DEBUGGING
                    if( idx >= m_count ) {
                        std::cerr << "Index out of range: " << idx << "; " << std::hex << first << "; " << second << std::dec << "; A size: " << N << "; B size: " << M << "; C size: " << c.m_size << std::endl;
                        assert( false );
                    }
#endif  // DEBUGGING

                    if( m_event_gen( idx ) ) {
                        sec |= bit_helper_type::bit_offset( idx );
                    }
                }

                block_type r = ((first & ~sec) | (second & sec));
                match_first = match_first && (r == first);
                match_second = match_second && (r == second );
                c_ptr[ i ] = r;

                j += bit_helper_type::BITS_PER_BLOCK;
                ++i;
            }

            if( match_first ) {
                return row_vector(a);
            } else if( match_second ) {
                return row_vector(b);
            } else {
                row_vector c( m_size, m_buffer );
                m_buffer = NULL;
                return c;
            }
        } else if( m_event_gen.getBaseSequence() != 0 ) {
            return row_vector( b );
        } else {
            return row_vector( a );
        }
    }

    virtual ~crossover_method() {
        if( m_buffer != NULL ) {
            delete [] m_buffer;
        }
    }

protected:
    event_generator_type m_event_gen;
    size_t m_size, m_count;

    block_type * m_buffer;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_CROSSOVER_METHOD_HPP_
