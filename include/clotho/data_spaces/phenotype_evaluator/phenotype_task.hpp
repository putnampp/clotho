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
#ifndef CLOTHO_PHENOTYPE_TASK_HPP_
#define CLOTHO_PHENOTYPE_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class AccumType >
class phenotype_task : public task {
public:

    typedef phenotype_task< BlockType, AccumType > self_type;

    typedef BlockType   block_type;
    typedef AccumType   accumulator_type;

    typedef clotho::utility::BitHelper< block_type >           bit_helper_type;
    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;

    phenotype_task( block_type * seq, unsigned int seq_len, accumulator_type & acc ) :
        m_seq( seq )
        , m_seq_length( seq_len )
        , m_accum( acc )
    {}

    phenotype_task( const self_type & other ) :
        m_seq( other.m_seq )
        , m_seq_length( other.m_seq_length )
        , m_accum( other.m_accum )
    {}

    void operator()() {

        size_t j = 0;
        for( unsigned int i = 0; i < m_seq_length; ++i ) {
            block_type b = m_seq[ i ];

            while( b ) {
                size_t b_idx = j + bit_walker_type::unset_next_index( b );
                m_accum( b_idx );
            }

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    virtual ~phenotype_task() {}
protected:
    block_type      * m_seq;
    unsigned int    m_seq_length;

    accumulator_type    m_accum;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PHENOTYPE_TASK_HPP_
