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
#ifndef CLOTHO_FREE_SPACE_MT_HPP_
#define CLOTHO_FREE_SPACE_MT_HPP_

#include <memory>

#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

#include "clotho/data_spaces/free_space/free_space_details.hpp"
#include "clotho/data_spaces/free_space/free_space_task.hpp"

#include "clotho/data_spaces/task/thread_pool.hpp"

namespace clotho {
namespace genetics {

template < class SizeType = unsigned int >
class FreeSpaceAnalyzerMT : public free_space_details< SizeType > {
public:
    typedef free_space_details< SizeType >  base_type;

    typedef SizeType            size_type;

    FreeSpaceAnalyzerMT( ) { }

    template < class SequenceSpaceType, class PoolType >
    void operator()( SequenceSpaceType & ss, PoolType & pool ) {
        typedef typename SequenceSpaceType::block_type block_type;

        this->resize( ss.column_count() );

        if ( ss.column_count() == 0 )  return;
        
        block_type * destF = new block_type[ 2 * ss.block_column_count() ];
        block_type * destV = destF + ss.block_column_count();

        memset( destF, 255,  ss.block_column_count() * sizeof(block_type) );
        memset( destV, 0,  ss.block_column_count() * sizeof(block_type) );

        process_space( ss.begin_row(0), destF, destV, ss.block_row_count(), ss.block_column_count(), pool );

        this->analyze_free_indices( destF, destV, ss.block_column_count(), ss.column_count() );

        delete [] destF;
    }

    virtual ~FreeSpaceAnalyzerMT() { }

protected:

    template < class BlockType, class PoolType >
    void process_space( BlockType * source, BlockType * destF, BlockType * destV, unsigned int block_rows, unsigned int block_columns, PoolType & pool ) {
        typedef free_space_task< BlockType > task_type;

        unsigned int tc = pool.pool_size() + 1; // + 1 for master thread

        // batching by columns allows this algorithm to elimnate the need for mutex locks
        unsigned int cpb = block_columns / tc;
        cpb += ((block_columns % tc > 0) ? 1 : 0);

        // cols will be decremented to cpb
        unsigned int cols = block_columns;
        while( cpb < cols ) {
            
            task_type t( source, destF, destV, block_rows, cpb, block_columns );

            pool.post( t );

            source += cpb;
            destF += cpb;
            destV += cpb;
            cols -= cpb;
        }

        if( 0 < cols ) {
            // last block being run on master thread
            task_type t(source, destF, destV, block_rows, cols, block_columns );
            t();
        }

        pool.sync();
    }
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class SizeType >
struct state_getter< clotho::genetics::FreeSpaceAnalyzerMT< SizeType > > {
    typedef clotho::genetics::FreeSpaceAnalyzerMT< SizeType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree fr;
        clotho::utility::add_value_array( fr, obj.free_begin(), obj.free_end() );
        
        s.put_child( "free", fr );
        s.add( "free_size", obj.free_size() );
        s.add( "variable_count", obj.variable_count() );
    }
};

}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_MT_HPP_
