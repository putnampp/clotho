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
#include "clotho/data_spaces/generators/buffered_small_crossover_event_generator.hpp"
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

    typedef small_crossover_event_generator< RNG, position_type >     event_generator_type;

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

template < class BlockType >
struct base_crossover_method {
    inline BlockType operator()( const BlockType & a, const BlockType & b, const BlockType & c) const {
        return ((a & ~c) | (b & c));
    }
};

template < class BlockType >
struct alt_crossover_method {
    inline BlockType operator()( const BlockType & a, const BlockType & b, const BlockType & c) const {
        return ((a & c) | (b & ~c));
    }
};

template < class RNG, class AlleleSpaceType, class BlockType >
class crossover_method< RNG, AlleleSpaceType, association_matrix< BlockType, row_grouped< 1 > > > {
public:
    typedef RNG                                                 random_engine_type;
    typedef AlleleSpaceType                                     allele_type;
    typedef association_matrix< BlockType, row_grouped< 1 > >   sequence_space_type;

    typedef typename sequence_space_type::raw_block_pointer     sequence_iterator;
    typedef typename sequence_space_type::raw_block_pointer     genome_iterator;

    typedef typename allele_type::position_type                 position_type;

    typedef buffered_small_crossover_event_generator< RNG, position_type >     event_generator_type;

    typedef typename sequence_space_type::block_type            block_type;
    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    typedef base_crossover_method< block_type >                 base_method;
    typedef alt_crossover_method< block_type >                  alt_method;

    crossover_method( RNG * rng, boost::property_tree::ptree & config ) :
        m_event_gen( rng, config )
    {}

    void initAlleles( allele_type & alls ) {
        m_event_gen.update( alls.getPositions(), alls.size() );
    }

/**
 *
 * @return soft max of result sequence
 */
    unsigned int crossover( sequence_iterator first, unsigned int N, sequence_iterator second, unsigned int M, sequence_iterator s ) {
        unsigned int res = 0;
        if( m_event_gen.generate() != 0 ) {
            if( m_event_gen.getBaseSequence() == 0 ) {
                base_method met;
                res = buffered_write_build_sequence( s, first, N, second, M, met );
            } else {
                alt_method met;
                res = buffered_write_build_sequence(s, first, N, second, M, met );
            }
        } else if( m_event_gen.getBaseSequence() != 0 ) {
            // use second strand
            res = copy_sequence( s, second, M );
        } else {
            res = copy_sequence( s, first, N );
        }

        return res;
    }

    virtual ~crossover_method() {}

protected:

/**
 *
 * Copies global sequence data to local stack buffers
 * Operates on stack buffers writing back to stack
 * Copies stack to final output vector
 *
 * Is read and write buffering beneficial?
 *
 * @return soft max of result sequence
 */
    template < class Method >
    unsigned int buffered_build_sequence( sequence_iterator res, sequence_iterator first, const unsigned int soft_N, sequence_iterator second, const unsigned int soft_M, const Method & m ) {
        const unsigned int N = 16;
        block_type first_buffer[ N ];   // 16 * 8 == 128 byte buffer
        block_type second_buffer[ N ];

        const unsigned int soft_max = (( soft_N > soft_M )? soft_N : soft_M );   // soft maximum length
        unsigned int last_block = 0;

        unsigned int l = 0;
        while( l < soft_max ) {
            unsigned int j = (soft_max - l);
            j = ((j > N) ? N : j);

            // read global data to local buffers
            memcpy( first_buffer, first + l, j * sizeof( block_type ) );
            memcpy( second_buffer, second + l, j * sizeof( block_type ) );

            unsigned int i = 0;
            while( i < j ) {
                const unsigned int offset = l + i;
                block_type first = ((offset < soft_N ) ? first_buffer[ i ] : bit_helper_type::ALL_UNSET);
                block_type second = ((offset < soft_M ) ? second_buffer[ i ] : bit_helper_type::ALL_UNSET);

                block_type mask = buffered_build_mask( first, second, offset * bit_helper_type::BITS_PER_BLOCK );

                mask = m( first, second, mask );
                second_buffer[ i ] = mask;
                ++i;
                if( mask ) {
                    last_block = offset + 1;
                }
            }

            // write to result sequence
            memcpy( res + l, second_buffer, j * sizeof(block_type) );
            l += j;
        }

        return last_block;
    }

/**
 * Unraveled the set bit decoding
 */
    block_type buffered_build_mask( block_type first, block_type second, unsigned int offset ) {
        block_type hets = first ^ second;
        block_type mask  = bit_helper_type::ALL_UNSET;   // mask state from second strand

        unsigned int het_buffer[ bit_helper_type::BITS_PER_BLOCK ]; // i.e. buffer_size = 64 * 4 = 256 bytes buffer

        if( hets ) {
            m_event_gen.update_buffer( offset );

            // decode set indices
            unsigned int j = 0;
            do {
                het_buffer[ j ] = bit_walker_type::unset_next_index( hets );
                ++j;
            } while( hets );

            // evaluate positions at relative indices
            unsigned int i = 0;
            while( i < j ) {
                unsigned int idx = het_buffer[ i ];
                if( m_event_gen( idx ) ) {
                    mask |= clotho::utility::bit_masks[ idx ];
                }
                ++i;
            }
        }

        return mask;
    }

/**
 *
 * Initial technique for building a crossover mask
 *
 *
 */
    block_type build_mask( block_type first, block_type second, unsigned int offset ) {
        block_type hets = first ^ second;
        block_type mask = bit_helper_type::ALL_UNSET;   // mask state from second strand

        if( hets ) {
            m_event_gen.update_buffer( offset );

            do {
                unsigned int idx = bit_walker_type::unset_next_index( hets );

                if( m_event_gen( idx ) ) {
                    mask |= clotho::utility::bit_masks[ idx ];
                }
            } while( hets );
        }

        return mask;
    }

    template < class Method >
    void build_sequence( sequence_iterator s, genome_iterator start, genome_iterator end, const Method & m ) {
        size_t j = 0;
        genome_iterator sec_ptr = end;
        while( start != end ) {
            block_type first = *start++;
            block_type second = *sec_ptr++;

            block_type mask = build_mask( first, second, j );

            *s++ = m(first, second, mask);

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

/**
 *
 * Uses soft max of both sequence
 *
 * Is write buffer beneficial?
 *
 * @return soft max of result sequence
 */
    template < class Method >
    unsigned int buffered_write_build_sequence( sequence_iterator res, sequence_iterator first, const unsigned int soft_N, sequence_iterator second, const unsigned int soft_M, const Method & m ) {

        const unsigned int N = 32;
        block_type  _buffer[ N ];   // 256 byte buffer

        const unsigned int soft_max = (( soft_N > soft_M )? soft_N : soft_M );   // soft maximum length
        unsigned int last_block = 0;

        unsigned int i = 0, j = 0;

        while( j < soft_max ) {
            block_type a = ((j < soft_N ) ? first[ j ] : bit_helper_type::ALL_UNSET);
            block_type b = ((j < soft_M ) ? second[ j ] : bit_helper_type::ALL_UNSET);

            block_type mask = build_mask( a, b, j * bit_helper_type::BITS_PER_BLOCK );

            mask = m( a, b, mask );
            _buffer[ i++ ] = mask;
            ++j;

            if( mask ) {
                last_block = j;
            }

            if( i >= N ) {
                // buffer is full
                memcpy( res, _buffer, N * sizeof(block_type) );
                res += N;
                i = 0;
            }
        }

        // write final buffer
        if( i ) {
            memcpy( res, _buffer, i * sizeof( block_type ) );
        }

        return last_block;
    }

/**
 *
 *
 * @return soft max of result sequence
 */
    template < class Method >
    unsigned int build_sequence( sequence_iterator res, sequence_iterator first, const unsigned int N, sequence_iterator second, const unsigned int M, const Method & m ) {

        const unsigned int W = (( N > M )? N : M );
        unsigned int i = 0;
        unsigned int last_block = 0;
        while( i < W ) {
            block_type a = ((i < N ) ? first[ i ] : bit_helper_type::ALL_UNSET);
            block_type b = ((i < M ) ? second[ i ] : bit_helper_type::ALL_UNSET);

            block_type mask = build_mask( a, b, i * bit_helper_type::BITS_PER_BLOCK );

            mask = m(a, b, mask );
            res[ i ] = mask;
            ++i;

            if( mask ) {
                last_block = i;
            }
        }

        return last_block;
    }

    unsigned int copy_sequence( sequence_iterator s, sequence_iterator start, const size_t W ) {
        memcpy( s, start, W * sizeof(block_type));
        return W;
    }

    event_generator_type m_event_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_CROSSOVER_METHOD_HPP_
