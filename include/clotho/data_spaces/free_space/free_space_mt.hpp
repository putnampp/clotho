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

namespace clotho {
namespace genetics {

template < class SizeType = unsigned int >
class FreeSpaceAnalyzerMT : public free_space_details< SizeType > {
public:

    typedef SizeType            size_type;

    FreeSpaceAnalyzerMT( ) { }

    template < class SequenceSpaceType >
    void operator()( SequenceSpaceType & ss, unsigned int tc ) {
        typedef typename SequenceSpaceType::block_type block_type;

        this->resize( ss.column_count() );
        
        block_type * destF = new block_type[ 2 * ss.block_column_count() ];
        block_type * destV = destF + ss.block_column_count();

        memset( destF, 0, 2 * ss.block_column_count() * sizeof(block_type) );

        process_space( ss.begin_row(0), destF, destV, ss.block_row_count(), ss.block_column_count(), tc );

        this->analyze_free_indices( destF, destV, ss.block_column_count(), ss.column_count() );

        delete [] destF;
    }

    virtual ~FreeSpaceAnalyzerMT() { }

protected:

    template < class BlockType >
    void process_space( BlockType * source, BlockType * destF, BlockType * destV, unsigned int block_rows, unsigned int block_columns, unsigned int tc ) {
        typedef free_space_task< BlockType > task_type;

        if( tc <= 1 ) {
            task_type task( source, destF, destV, block_rows, block_columns );
            task();
        } else {
            boost::asio::io_service service;
            boost::thread_group threads;
            boost::mutex mutef, mutev;

            {
                std::unique_ptr< boost::asio::io_service::work > work( new boost::asio::io_service::work( service ) );
                if( tc > block_rows ) {
                    tc = block_rows;
                }

                for( size_type i = 0; i < tc; ++i )
                    threads.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );

                unsigned int rpb = block_rows / tc;
                rpb += ((block_rows % tc > 0) ? 1 : 0);

                while( block_rows ) {
                    unsigned int rows = (( rpb < block_rows ) ? rpb : block_rows );
                    std::shared_ptr< task_type > t( new task_type( source, destF, destV, rows, block_columns ) );
                    
                    service.post( boost::bind( &task_type::evaluate_ts, t, &mutef, &mutev) );
                    block_rows -= rows;
                    source += rows * block_columns;
                }
            }

            threads.join_all();
        }
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
