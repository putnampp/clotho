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
#ifndef CLOTHO_COLUMN_ALIGNED_CROSSOVER_METHOD_HPP_
#define CLOTHO_COLUMN_ALIGNED_CROSSOVER_METHOD_HPP_

#include "clotho/data_spaces/crossover/crossover_method_def.hpp"
#include "clotho/data_spaces/association_matrix/column_aligned_association_matrix.hpp"

#include "clotho/data_spaces/generators/crossover_event_generator.hpp"
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class AlleleSpaceType, class BlockType >
class crossover_method< RNG, AlleleSpaceType, association_matrix< BlockType, column_aligned > > {
public:
    typedef RNG                                                 random_engine_type;
    typedef AlleleSpaceType                                     allele_type;
    typedef association_matrix< BlockType, column_aligned >     sequence_space_type;

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
        m_event_gen.update( alls.position_begin(), alls.position_end() );
    }

    void crossover( genome_iterator start, genome_iterator end, sequence_iterator s, size_t p_step, size_t s_step ) {

        size_t N = m_event_gen.generate();

        if( N != 0 ) {
            p_step--;
            size_t j = 0;
            while( start != end ) {
                block_type first = *start++;
                block_type second = *start;

                block_type hets = first ^ second;
                block_type sec  = bit_helper_type::ALL_UNSET;   // mask state from second strand

                while( hets ) {
                    size_t idx = bit_walker_type::unset_next_index( hets ) + j;

                    if( m_event_gen( idx ) ) {
                        sec |= bit_helper_type::bit_offset( idx );
                    }
                }

                hets = ((first & ~sec) | (second & sec));
                *s = hets;

                j += bit_helper_type::BITS_PER_BLOCK;

                start += p_step;
                s += s_step;
            }
        } else {
            if( m_event_gen.getBaseSequence() != 0 ) {
                ++start;
                ++end;
            }

            while( start != end ) {
                *s = *start;

                start += p_step;
                s += s_step;
            }
        }
    }
    virtual ~crossover_method() {}

protected:
    event_generator_type m_event_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_COLUMN_ALIGNED_CROSSOVER_METHOD_HPP_
