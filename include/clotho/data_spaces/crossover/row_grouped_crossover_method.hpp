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
        m_event_gen.update( alls.getPositions(), alls.size() );
    }

    void crossover( genome_iterator start, genome_iterator end, sequence_iterator s, size_t p_step, size_t s_step ) {

        size_t N = m_event_gen.generate();

        genome_iterator sec_ptr = end;
        if( N != 0 ) {
            size_t j = 0;
            while( start != end ) {
                block_type first = *start++;
                block_type second = *sec_ptr++;

                block_type hets = first ^ second;
                block_type sec  = bit_helper_type::ALL_UNSET;   // mask state from second strand

                while( hets ) {
                    size_t idx = bit_walker_type::unset_next_index( hets ) + j;

                    if( m_event_gen( idx ) ) {
                        sec |= bit_helper_type::bit_offset( idx );
                    }
                }

                *s++ = ((first & ~sec) | (second & sec));

                j += bit_helper_type::BITS_PER_BLOCK;
            }
        } else {
            size_t N = (end - start);

            if( m_event_gen.getBaseSequence() != 0 ) {
                // use second strand
                start = end;
            }

            for( size_t i = 0; i < N; ++i )
                s[i] = start[i];
//            size_t i = 0;
//            size_t M = N % 4;
//            while( i < M ) {
//                s[i] = start[i];
//                ++i;
//            }
//
//            while( i < N ) {
//                s[ i ] = start[ i ];
//                s[ i + 1 ] = start[i + 1];
//                s[ i + 2 ] = start[i + 2 ];
//                s[ i + 3 ] = start[i + 3 ];
//                i += 4;
//            }
//
        }
    }

    virtual ~crossover_method() {}

protected:
    event_generator_type m_event_gen;
};
}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_GROUPED_CROSSOVER_METHOD_HPP_
