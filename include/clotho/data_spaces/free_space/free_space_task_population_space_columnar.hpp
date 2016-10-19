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
#ifndef CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_COLUMNAR_HPP_
#define CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_COLUMNAR_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/free_space/free_space_task.hpp"
#include "clotho/data_spaces/population_space/population_space_columnar.hpp"

namespace clotho {
namespace genetics {


template < class BlockType, class WeightType >
class free_space_task < population_space_columnar< BlockType, WeightType > > : public task {
public:
    typedef free_space_task< population_space_columnar< BlockType, WeightType > > self_type;

    typedef population_space_columnar< BlockType, WeightType >       space_type;
    typedef typename space_type::block_type                          block_type;

    typedef typename space_type::row_pointer                         row_pointer;

    typedef typename space_type::bit_helper_type                     bit_helper_type;

    free_space_task( space_type * ss, block_type * fixed_dest, block_type * var_dest, unsigned int block_start, unsigned int block_end ) :
        m_space(ss)
        , m_destF( fixed_dest + block_start )
        , m_destF_end( fixed_dest + block_end )
        , m_destV( var_dest + block_start )
        , m_start( block_start )
        , m_end( block_end )
    {}

    free_space_task( const self_type & other ) :
        m_space( other.m_space )
        , m_destF( other.m_destF )
        , m_destF_end( other.m_destF_end )
        , m_destV( other.m_destV )
        , m_start( other.m_start )
        , m_end( other.m_end )
    {}

    void operator()() {

        unsigned int i = m_start;
        while( i < m_end ) {
            block_type fx = bit_helper_type::ALL_SET;
            block_type var = bit_helper_type::ALL_UNSET;

            row_pointer first = m_space->begin_block_row( i ), last = m_space->end_block_row(i);
            while( first != last ) {
                fx &= *first;
                var |= *first;

                ++first;
            }

            m_destF[ i ] = fx;
            m_destV[ i ] = var;

            ++i;
        }

    }

    virtual ~free_space_task() {}

protected:
    space_type  * m_space;

    block_type  * m_destF, * m_destF_end, * m_destV;
    unsigned int m_start, m_end;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_COLUMNAR_HPP_
