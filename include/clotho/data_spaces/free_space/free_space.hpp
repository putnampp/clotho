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
#ifndef CLOTHO_FREE_SPACE_HPP_
#define CLOTHO_FREE_SPACE_HPP_

#include <vector>
#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
class FreeSpaceAnalyzer {
public:
    typedef GeneticSpaceType    genetic_space_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef typename genetic_space_type::block_iterator block_iterator;

    typedef typename association_type::block_type   block_type;

    typedef std::vector< size_t >                   index_vector;
    typedef typename free_index_vector::iterator    iterator;
    typedef typename free_index_vector::const_iterator    const_iterator;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;

    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    FreeSpaceAnalyzer(){}

    FreeSpaceAnalyzer( genetic_space_type & gs ) {
        update( gs );
    }

    void update( genetic_space_type & gs ) {
        m_free.clear();
        m_fixed.clear();
        m_lost.clear();

        typedef std::vector< block_type >   block_vector;

        block_iterator b_it = gs.getBlockIterator();

        size_t M = gs.allele_count();

        block_type fx = bit_helper_type::ALL_SET, var = bit_helper_type::ALL_UNSET;
        size_t i = 0, j = 0;
        while( b_it.hasNext() ) {
            block_type v = b.next();
            
            fx &= v;
            var |= v;

            if( ++i >= gs.block_column_count() ) {
                block_type ls = ~(fx | var);

                while( fx ) {
                    size_t idx = bit_walker_type::unset_next_index( fx ) + j;
                    if( idx < M ) {
                        m_fixed.push_back( idx);
                        m_free.push_back( idx );
                    }
                }

                while( ls ) {
                    size_t idx = bit_walker_type::unset_next_index( ls ) + j;
                    if( idx < M ) {
                        m_lost.push_back( idx );
                        m_free.push_back( idx );
                    }
                }

                fx = bit_helper_type::ALL_SET;
                var = bit_helper_type::ALL_UNSET;

                j += bit_helper_type::BITS_PER_BLOCK;
            }
        }

        assert( i == gs.block_row_count() );
        assert( _fixed.size() == i );
    }

    size_t free_size() const {
        return m_free.size();
    }

    iterator free_begin() {
        return m_free.begin();
    }

    iterator free_end() {
        return m_free.end();
    }

    const_iterator free_begin() const {
        return m_free.begin();
    }

    const_iterator free_end() const {
        return m_free.end();
    }

    size_t fixed_size() const {
        return m_fixed.size();
    }

    iterator fixed_begin() {
        return m_fixed.begin();
    }

    iterator fixed_end() {
        return m_fixed.end();
    }

    const_iterator fixed_begin() const {
        return m_fixed.begin();
    }

    const_iterator fixed_end() const {
        return m_fixed.end();
    }

    size_t lost_size() const {
        return m_lost.size();
    }

    iterator lost_begin() {
        return m_lost.begin();
    }

    iterator lost_end() {
        return m_lost.end();
    }

    const_iterator lost_begin() const {
        return m_lost.begin();
    }

    const_iterator lost_end() const {
        return m_lost.end();
    }

    virtual ~FreeSpaceAnalyzer() {}
protected:
    index_vector m_free, m_fixed, m_lost;

};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_HPP_
