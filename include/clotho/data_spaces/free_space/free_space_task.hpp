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
#ifndef CLOTHO_FREE_SPACE_TASK_HPP_
#define CLOTHO_FREE_SPACE_TASK_HPP_

#include "clotho/data_spaces/task/task.hpp"

#include <boost/thread/thread.hpp>

namespace clotho {
namespace genetics {

template < class BlockType >
class free_space_task : public task {
public:
    typedef BlockType   block_type;

    typedef free_space_task< BlockType > self_type;

    free_space_task( block_type * src, block_type * fixed_dest, block_type * var_dest, unsigned int block_rows, unsigned int block_cols, unsigned int row_width ) :
        m_source( src )
        , m_fixed_dest( fixed_dest )
        , m_variable_dest( var_dest )
        , m_block_rows( block_rows )
        , m_block_columns( block_cols )
        , m_row_width( row_width )
    {}


    free_space_task( const self_type & other ) :
        m_source( other.m_source )
        , m_fixed_dest( other.m_fixed_dest )
        , m_variable_dest( other.m_variable_dest )
        , m_block_rows( other.m_block_rows )
        , m_block_columns( other.m_block_columns )
        , m_row_width( other.m_row_width )
    {}

    void evaluate_ts( boost::mutex * flock, boost::mutex * vlock) {

        block_type * tempF = new block_type[ 2 * m_block_columns ];
        block_type * tempV = tempF + m_block_columns;

        memset( tempF, 0, 2 * m_block_columns * sizeof( block_type ) );

        evaluate( tempF, tempV );

        {
            // lock destination
            boost::lock_guard< boost::mutex > l( *flock );

            for( unsigned int i = 0; i < m_block_columns; ++i ) {
                m_fixed_dest[ i ] &= tempF[ i ];
            }
        }

        {
            // lock destination
            boost::lock_guard< boost::mutex > l( *vlock );
            for( unsigned int i = 0; i < m_block_columns; ++i ) {
                m_variable_dest[ i ] |= tempV[ i ];
            }
        }

        delete [] tempF;
    }

    void operator()( ) {
        evaluate( m_fixed_dest, m_variable_dest );
    }

    void evaluate( block_type * tempF, block_type * tempV ) {
#ifdef DEBUGGING
        BOOST_LOG_TRIVIAL(debug) << "Evaluating Free Space: " << m_block_rows << " x " << m_block_columns;
#endif  // DEBUGGING

        block_type * src = m_source;
        unsigned int offset = 0;
        for( unsigned int i = 0; i < m_block_rows; ++i ) {
            for( unsigned int j = 0; j < m_block_columns; ++j ) {
                block_type b = src[ offset + j ];
                tempF[ j ] &= b;
                tempV[ j ] |= b;
            }
            offset += m_row_width;
        }
    }

    virtual ~free_space_task() {}

protected:
    block_type * m_source, * m_fixed_dest, * m_variable_dest;
    unsigned int m_block_rows, m_block_columns, m_row_width;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_TASK_HPP_

