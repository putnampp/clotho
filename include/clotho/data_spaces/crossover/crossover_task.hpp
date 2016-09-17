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
#ifndef CLOTHO_CROSSOVER_TASK_HPP_
#define CLOTHO_CROSSOVER_TASK_HPP_

#include "clotho/data_spaces/generators/small_crossover_event_generator.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

#include "clotho/data_spaces/task/task.hpp"
#include <iostream>

#include <boost/thread/thread.hpp>
#include <boost/property_tree/ptree.hpp>

#include "clotho/data_spaces/crossover/crossover_details.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PositionType, class BlockType >
class crossover_task : public task {
public:
    typedef crossover_task< RNG, PositionType, BlockType > self_type;

    typedef RNG         random_generator_type;
    typedef BlockType   block_type;

    typedef small_crossover_event_generator< RNG, PositionType > event_generator_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    typedef base_crossover_method< block_type > base_method;
    typedef alt_crossover_method< block_type > alt_method;

    template < typename SeedSeq >
    crossover_task( SeedSeq & seed
        , double rho
        , double bias
        , PositionType * pos
        , unsigned int Npos
        , block_type * p0
        , unsigned int p0_soft_max
        , block_type * p1
        , unsigned int p1_soft_max
        , block_type * offspring
        , unsigned int * offspring_soft_max
        ) :
        m_rng( seed )
        , m_event_gen( &m_rng, rho, bias, pos, Npos )
        , m_p0( p0 )
        , m_p0_soft_max( p0_soft_max )
        , m_p1( p1 )
        , m_p1_soft_max( p1_soft_max )
        , m_offspring( offspring )
        , m_offspring_soft_max( offspring_soft_max )
    {}

    crossover_task( const self_type & other ) :
        m_rng( other.m_rng )
        , m_event_gen( other.m_event_gen )
        , m_p0( other.m_p0 )
        , m_p0_soft_max( other.m_p0_soft_max )
        , m_p1( other.m_p1 )
        , m_p1_soft_max( other.m_p1_soft_max )
        , m_offspring( other.m_offspring )
        , m_offspring_soft_max( other.m_offspring_soft_max )
    {}

    void operator()() {
        size_t N = m_event_gen.generate();

        if( N != 0 ) {
            if( m_event_gen.getBaseSequence() == 0 ) {
                base_method met;
                build_sequence( met );
            } else {
                alt_method met;
                build_sequence( met );
            }
        } else if( m_event_gen.getBaseSequence() == 0 ) {
            copy_sequence( m_p0, m_p0_soft_max );
        } else {
            copy_sequence( m_p1, m_p1_soft_max );
        }
    }

    virtual ~crossover_task() {}

protected:

    template < class Method >
    void build_sequence( const Method & m ) {
        const unsigned int W = (( m_p0_soft_max < m_p1_soft_max )? m_p0_soft_max : m_p1_soft_max );
        unsigned int i = 0;
        unsigned int last_block = 0;

        while( i < W ) {
            block_type a = m_p0[ i ];
            block_type b = m_p1[ i ];

            block_type mask = build_mask( a, b, i * bit_helper_type::BITS_PER_BLOCK );

            mask = m(a, b, mask );
            m_offspring[ i ] = mask;
            ++i;

            if( mask ) {
                last_block = i;
            }
        }

        while( i < m_p0_soft_max ) {
            block_type a = m_p0[i];

            block_type mask = build_mask( a, bit_helper_type::ALL_UNSET, i * bit_helper_type::BITS_PER_BLOCK );

            mask = m(a, bit_helper_type::ALL_UNSET, mask );
            m_offspring[ i ] = mask;
            ++i;

            if( mask ) {
                last_block = i;
            }
        }

        while( i < m_p1_soft_max ) {
            block_type b = m_p1[i];

            block_type mask = build_mask( bit_helper_type::ALL_UNSET, b, i * bit_helper_type::BITS_PER_BLOCK );

            mask = m( bit_helper_type::ALL_UNSET, b, mask );
            m_offspring[ i ] = mask;
            ++i;

            if( mask ) {
                last_block = i;
            }
        }

        *m_offspring_soft_max = last_block;
    }

    inline void copy_sequence( block_type * src, const size_t W ) {
        memcpy( m_offspring, src, W * sizeof( block_type ) );
        *m_offspring_soft_max = W;
    }

    block_type build_mask( block_type first, block_type second, unsigned int offset ) {
        block_type hets = first ^ second;
        block_type mask = bit_helper_type::ALL_UNSET;   // mask state from m_p1 strand

        if( hets ) {
            do {
                unsigned int idx = bit_walker_type::unset_next_index( hets );

                if( m_event_gen( offset + idx ) ) {
                    mask |= clotho::utility::bit_masks[ idx ];
                }
            } while( hets );
        }

        return mask;
    }

    random_generator_type m_rng;
    event_generator_type m_event_gen;

    block_type * m_p0;
    unsigned int m_p0_soft_max;

    block_type * m_p1;
    unsigned int m_p1_soft_max;

    block_type * m_offspring;
    unsigned int * m_offspring_soft_max;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_TASK_HPP_
