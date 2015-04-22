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
#ifndef FITNESS_BITSET_HPP_
#define FITNESS_BITSET_HPP_

#include <boost/dynamic_bitset.hpp>
#include "utility/lowest_bit.h"

#include "accessor.hpp"

namespace clotho {
namespace fitness {

template < class Block, class Allocator, class Alphabet, class HomPolicy, class HetPolicy, class ResultType >
class fitness_bitset;

namespace boost {

template < class Block, class Allocator, class Alphabet, class HomPolicy, class HetPolicy, class ResultType >
inline void to_block_range( dynamic_bitset< Block, Allocator > & alt, fitness_bitset< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType > fitness ) {
    typedef typename std::vector< Block, Allocator >::iterator iterator;
    typedef Block block_type;

    iterator first = fitness.m_base->m_bits.begin(), last = fitness.m_base->m_bits.end(),
             _first = alt.m_bits.begin(), _last = alt.m_bits.end();

    fitness( first, last, _first, _last );
}

}   // namespace boost

template < class Block, class Allocator, class Alphabet, class HomPolicy, class HetPolicy, class ResultType >
class fitness_bitset {
public:
    typedef fitness_bitset< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType > self_type;

    typedef boost::dynamic_bitset< Block, Allocator >   bitset_type;
    typedef ResultType                                  result_type;

    typedef typename Alphabet::active_iterator          active_iterator;
    typedef lowest_bit_64K                              lowest_bit_map;

    friend void boost::to_block_range< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType >( boost::dynamic_bitset< Block, Allocator > &, fitness_bitset< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType > );

    fitness_bitset( bitset_type * base, typename Alphabet::pointer alpha, HomPolicy & homozygous, HetPolicy & heterozygous, result_type initial ) :
        m_base( base ), m_alpha( alpha ), m_hom(&homozygous), m_het(&heterozygous), m_result( initial) {
    }

    fitness_bitset( const self_type & other ) :
        m_base( other.m_base ), m_alpha( other.m_alpha ), m_hom(other.m_hom), m_het( other.m_het ), m_result( other.m_result ) {
    }

    template < class BlockIterator >
    void operator()( BlockIterator base_first, BlockIterator base_last, BlockIterator alt_first, BlockIterator alt_last ) {
        active_iterator seq_pos = m_alpha->active_begin();

        while( true ) {
            if( base_first == base_last ) {
                while( alt_first != alt_last ) {
                    bit_walker_unrolled( (*alt_first++), seq_pos, m_het);
                    seq_pos += bitset_type::bits_per_block;
                }
                break;
            }

            if( alt_first == alt_last ) {
                while( base_first != base_last ) {
                    bit_walker_unrolled( (*base_first++), seq_pos, m_het );
                    seq_pos += bitset_type::bits_per_block;
                }
                break;
            }
            Block base = (*base_first++), alt = (*alt_first++);

            bit_walker_unrolled( base & alt, seq_pos, m_hom );
            bit_walker_unrolled( base ^ alt, seq_pos, m_het );

            seq_pos += bitset_type::bits_per_block;
        }
    }

    bitset_type * getBase() {
        return m_base;
    }

    result_type getResult() {
        return m_result;
    }

    virtual ~fitness_bitset() {}
protected:

    inline void walk_block_bits( Block base, Block alt, unsigned int pos_offset ) {
        Block bits = (base | alt);      // union
        Block hets = (base ^ alt );

        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            // NOTE: It may be more efficient to subdivide a block into byte regions
            // to allow for skipping of long regions of 0. Basically, recursively
            // walk the set bits by regions
            if( bits & 1 ) {    // if the locus is heterozygous
                // determine from which recombination region this locus is inherited from (base or alternate)
                if( hets & 1 ) {
                    m_het( m_result, *(*m_alpha)[ pos_offset ]->second.first);
                } else {
                    m_hom( m_result, *(*m_alpha)[ pos_offset ]->second.first );
                }
            }
            bits >>= 1;
            hets >>= 1;
            ++pos_offset;
        }
    }

    inline void walk_block_bits_sparse( Block base, Block alt, unsigned int pos_offset ) {
        // performance should degrade as span of 0's is less than or equal to 8
        Block bits = (base | alt);      // symmetric difference
        Block hets = (base ^ alt );

        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            if( bits & 255 ) {
                unsigned int n = 8;
                while( n-- ) {
                    if( bits & 1 ) {    // if the locus is heterozygous
                        // determine from which recombination region this locus is inherited from (base or alternate)
                        if( hets & 1 ) {
                            m_het( m_result, *(*m_alpha)[ pos_offset ]->second.first);
                        } else {
                            m_hom( m_result, *(*m_alpha)[ pos_offset ]->second.first );
                        }
                        n = 0;  // break out of loop to recheck state
                    }
                    bits >>= 1;
                    hets >>= 1;
                    ++pos_offset;
                }
            } else {
                bits >>= 8;
                hets >>= 8;
                pos_offset += 8;
            }
        }
    }

    inline void walk_block_bits_sparse2( unsigned long base, unsigned long alt, unsigned int pos_offset ) {
        // performance should degrade as span of 0's is less than or equal to 8
        Block bits = (base | alt);      // symmetric difference
        Block hets = (base ^ alt);

        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            if( bits & 0x00000000000000FF ) {
                unsigned int n = 8;
                while( n-- ) {
                    if( bits & 1 ) {    // if the locus is heterozygous
                        // determine from which recombination region this locus is inherited from (base or alternate)
                        if( hets & 1 ) {
                            m_het( m_result, *(*m_alpha)[ pos_offset ]->second.first);
                        } else {
                            m_hom( m_result, *(*m_alpha)[ pos_offset ]->second.first );
                        }
                        n = 0;  // break out of loop to recheck state
                    }
                    bits >>= 1;
                    hets >>= 1;
                    ++pos_offset;
                }
            } else if( bits & 0x000000000000FF00) {
                bits >>= 8;
                hets >>= 8;
                pos_offset += 8;
            } else if( bits & 0x00000000FFFF0000) {
                bits >>= 16;
                hets >>= 16;
                pos_offset += 16;
            } else {
                bits >>= 32;
                hets >>= 32;
                pos_offset += 32;
            }
        }
    }

    inline void walk_block_bits_sparse3( unsigned long base, unsigned long alt, unsigned int pos_offset ) {
        bit_walker( base & alt, pos_offset, &m_hom );
        bit_walker( base ^ alt, pos_offset, &m_het );
    }

    /*
        template < class OP >
        inline void bit_walker( Block _bits, unsigned int pos_offset, OP * op ) {
    //        unsigned int res_offset = 0;
            while( _bits ) {
                unsigned char low_byte = (unsigned char)(_bits & 0x00000000000000FF);
    //            unsigned int _offset = res_offset;
                unsigned int _offset = pos_offset;
                while( low_byte ) {
                    const lowest_bit_map::value_type & v = low_bit_map[low_byte];

                    (*op)( m_result, *(*m_alpha)[ _offset + v.bit_index ]->second.first );

                    low_byte = v.next;
                    _offset += v.bit_shift_next;
                }

                _bits >>= 8;
                pos_offset += 8;
    //            res_offset += 8;
            }
        }*/

    template < class OP >
    inline void bit_walk( typename lowest_bit_map::block_type low_byte, active_iterator base_it, OP * op, unsigned int _offset ) {
        const lowest_bit_map::value_type * v = low_bit_map.begin() + low_byte;
        do {
            (*op)( m_result, accessor::get< typename Alphabet::active_iterator, typename Alphabet::allele_t >(base_it + _offset + v->bit_index));

            _offset += v->bit_shift_next;
            v = v->next_ptr;
        } while( v != NULL );
    }

    template < class OP >
    inline void block_walker( unsigned long _bits, active_iterator base_it, OP * op, unsigned short t ) {
        typename lowest_bit_map::block_type low_byte = (typename lowest_bit_map::block_type)(_bits);

        if( low_byte ) {    // block 0
            bit_walk( low_byte, base_it, op, 0 );
        }
        _bits /= lowest_bit_map::max_values;

        if( !_bits ) return;

        low_byte = (typename lowest_bit_map::block_type)(_bits);

        if( low_byte ) {    // block 1
            bit_walk( low_byte, base_it, op, lowest_bit_map::block_width );
        }
        _bits /= lowest_bit_map::max_values;

        if( !_bits ) return;
        low_byte = (typename lowest_bit_map::block_type)(_bits);

        if( low_byte ) {    // block 2
            bit_walk( low_byte, base_it, op, 2 * lowest_bit_map::block_width );
        }
        _bits /= lowest_bit_map::max_values;

        low_byte = (typename lowest_bit_map::block_type)(_bits);

        if( low_byte ) {    // block 3
            bit_walk( low_byte, base_it, op, 3 * lowest_bit_map::block_width );
        }
    }

    template < class OP >
    inline void bit_walker_unrolled( Block _bits, active_iterator base_it, OP * op ) {
        if( !_bits ) return;

        block_walker( _bits, base_it, op, (typename lowest_bit_map::block_type) 0 );
    }

    template < class OP >
    inline void bit_walker( Block _bits, active_iterator base_it , OP * op ) {
        unsigned int res_offset = 0;
        while( _bits ) {
            typename lowest_bit_map::block_type low_byte = (typename lowest_bit_map::block_type)(_bits);
            if( low_byte ) {
                unsigned int _offset = res_offset;
                const lowest_bit_map::value_type * v = low_bit_map.begin() + low_byte;
                do {
                    //(*op)( m_result, *(*(base_it + _offset + v->bit_index))->second.first );
                    (*op)( m_result, accessor::get< typename Alphabet::active_iterator, typename Alphabet::allele_t >(base_it + _offset + v->bit_index));

                    _offset += v->bit_shift_next;
                    v = v->next_ptr;
                } while( v != NULL );
            }

            _bits >>= lowest_bit_map::block_width;
            res_offset += lowest_bit_map::block_width;
        }
    }

    template < class OP >
    inline void rec_bit_walker( unsigned short bits, active_iterator base_it, OP * op, unsigned int res_offset ) {
        if( bits == 0 ) return;

        const lowest_bit_map::value_type * v = low_bit_map.begin() + bits;
        do {
            (*op)( m_result, *(*(base_it + res_offset + v->bit_index))->second.first );

            res_offset += v->bit_shift_next;
            v = v->next_ptr;
        } while( v != NULL );
    }

    template < class OP >
    inline void rec_bit_walker( unsigned int bits, active_iterator pos_offset, OP * op, unsigned int res_offset ) {
        if( bits == 0 ) return;

        static const unsigned char WIDTH = 16;

        unsigned short _bits = (unsigned short) bits;
        rec_bit_walker( _bits, pos_offset, op, res_offset );

        _bits = (unsigned short)(bits >> WIDTH);

        rec_bit_walker( _bits, pos_offset, op, res_offset + WIDTH );
    }

    template < class OP >
    inline void rec_bit_walker( unsigned long bits, active_iterator pos_offset, OP * op, unsigned int res_offset = 0 ) {
        if( bits == 0 ) return;

        static const unsigned char WIDTH = 32;

        unsigned int _bits = (unsigned int) bits;
        rec_bit_walker( _bits, pos_offset, op, res_offset );

        _bits = (unsigned int)(bits >> WIDTH);

        rec_bit_walker( _bits, pos_offset, op, res_offset + WIDTH );
    }

    bitset_type * m_base;
    typename Alphabet::pointer  m_alpha;
    HomPolicy * m_hom;
    HetPolicy * m_het;
    result_type m_result;

    static const lowest_bit_map low_bit_map;
};


template < class Block, class Allocator, class Alphabet, class HomPolicy, class HetPolicy, class ResultType >
const typename fitness_bitset< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType >::lowest_bit_map
fitness_bitset< Block, Allocator, Alphabet, HomPolicy, HetPolicy, ResultType >::low_bit_map;

}   // namespace fitness
}   // namespace clotho

#endif  // FITNESS_BITSET_HPP_
