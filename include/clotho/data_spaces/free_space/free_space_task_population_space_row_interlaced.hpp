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
#ifndef CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_ROW_HPP_
#define CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/free_space/free_space_task.hpp"
#include "clotho/data_spaces/population_space/population_space_row.hpp"

namespace clotho {
namespace genetics {


template < class BlockType, class WeightType >
class free_space_task < population_space_row< BlockType, WeightType > > : public task {
public:
    typedef free_space_task< population_space_row< BlockType, WeightType > > self_type;

    typedef population_space_row< BlockType, WeightType >       space_type;
    typedef typename space_type::block_type                          block_type;

    typedef typename space_type::genome_pointer                      genome_pointer;

    typedef typename space_type::bit_helper_type                     bit_helper_type;

    free_space_task( space_type * ss, block_type * fixed_dest, unsigned int block_start, unsigned int block_end ) :
        m_space(ss)
        , m_destF( fixed_dest + 2 * block_start )
        , m_start( block_start )
        , m_end( block_end )
    {}

    free_space_task( const self_type & other ) :
        m_space( other.m_space )
        , m_destF( other.m_destF )
        , m_start( other.m_start )
        , m_end( other.m_end )
    {}

    void operator()() {
        assert( m_start < m_end );

        const unsigned int N = m_space->haploid_genome_count();
        //const unsigned int M = (( m_end < m_space->getMaxBlocks() ) ? m_end : m_space->getMaxBlocks());

        for( unsigned i = 1; i < N; ++i ) {
            genome_pointer first = m_space->begin_genome( i );
            genome_pointer last = first + m_end;
            first += m_start;

            // m_destF and m_destV have already been shifted according to the starting block
            block_type * df = m_destF;

            while( first != last ) {
                *df++ &= *first;
                *df++ |= *first;
                ++first;
            }
        }

    }

    virtual ~free_space_task() {}

protected:
    space_type  * m_space;

    block_type  * m_destF;
    unsigned int m_start, m_end;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_ROW_HPP_
