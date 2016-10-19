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
#ifndef CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_COLUMNAR_HPP_
#define CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_COLUMNAR_HPP_

#include "clotho/data_spaces/free_space/free_space_mt.hpp"
#include "clotho/data_spaces/population_space/population_space_columnar.hpp"
#include "clotho/data_spaces/free_space/free_space_tasks.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
class FreeSpaceAnalyzerMT < population_space_columnar< BlockType, WeightType >, SizeType > : public free_space_details< SizeType > {
public:

    typedef free_space_details< SizeType >              base_type;

    typedef population_space_columnar< BlockType, WeightType >   space_type;
    typedef typename space_type::block_type block_type;

    typedef free_space_task< space_type >    task_type;

    FreeSpaceAnalyzerMT() {}

    template < class PoolType >
    void operator()( space_type * ss, PoolType & pool ) {
        
        const unsigned int A = ss->getMaxAlleles();
        const unsigned int B = ss->getMaxBlocks();

        if( A == 0 ) return;

        this->resize( A );

        block_type * destF = new block_type[ 2 * B ];
        block_type * destV = destF + B;

        process_space( ss, destF, destV, B, pool );

        this->analyze_free_indices( destF, destV, B, A );

        delete [] destF;
    }

    virtual ~FreeSpaceAnalyzerMT() {}

protected:

    template < class PoolType >
    void process_space( space_type * source, block_type * destF, block_type * destV, unsigned int block_columns, PoolType & pool ) {

        const unsigned int tc = pool.pool_size() + 1; // + 1 for master thread

        // batching by columns allows this algorithm to elimnate the need for mutex locks
        const unsigned int cpb = block_columns / tc + ((block_columns % tc > 0) ? 1 : 0);

        // cols will be decremented to cpb
        // unsigned int cols = block_columns;
        unsigned int cols = 0;
        while( cols + cpb < block_columns  ) {
            
            task_type t( source, destF, destV, cols, cols + cpb );

            pool.post( t );
            cols += cpb;
        }

        if( cols < block_columns ) {
            // last block being run on master thread
            task_type t(source, destF, destV,  cols, block_columns );
            t();
        }

        pool.sync();
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_COLUMNAR_HPP_

