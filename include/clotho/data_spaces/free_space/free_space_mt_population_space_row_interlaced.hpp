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
#ifndef CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_ROW_HPP_
#define CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/free_space/free_space_mt.hpp"
#include "clotho/data_spaces/population_space/population_space_row.hpp"
#include "clotho/data_spaces/free_space/free_space_tasks.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType, class SizeType >
class FreeSpaceAnalyzerMT < population_space_row< BlockType, WeightType >, SizeType > : public free_space_details< SizeType > {
public:

    typedef free_space_details< SizeType >              base_type;

    typedef population_space_row< BlockType, WeightType >   space_type;
    typedef typename space_type::block_type                 block_type;
    typedef typename space_type::bit_helper_type            bit_helper_type;

    typedef typename space_type::genome_pointer             genome_pointer;

    typedef free_space_task< space_type >                   task_type;

    FreeSpaceAnalyzerMT() {}

    template < class PoolType >
    void operator()( space_type * ss, PoolType & pool ) {
        
        const unsigned int A = ss->getMaxAlleles();
        const unsigned int B = ss->getMaxBlocks();

        if( A == 0 ) return;

        this->resize( A );

        block_type * destF = new block_type[ 2 * B ];

        if( ss->haploid_genome_count() == 0 ) {
            // there are no genomes; so there are no fixed or variable alleles
            // and no need to process the data space either
            memset( destF, 0, 2 * B * sizeof(block_type));
        } else {
            // copy the first sequence as the start of the free space logic
            genome_pointer first = ss->begin_genome(0), last = ss->end_genome(0);
            assert( (last - first) == B );

            block_type * df = destF;
            while( first != last ) {
                *df++ = *first;
                *df++ = *first;
                ++first;
            }

            if( (A / bit_helper_type::BITS_PER_BLOCK) != B) {
                block_type * s = destF +  2 * (A / bit_helper_type::BITS_PER_BLOCK);
                block_type * e = destF + 2 * B;
                // account for padding in tail blocks
                block_type padding_mask = (bit_helper_type::ALL_SET << (A % bit_helper_type::BITS_PER_BLOCK));

                *s++ &= ~padding_mask;
                *s++ |= padding_mask;
                
                while( s != e ) {
                    *s++ = bit_helper_type::ALL_UNSET;
                    *s++ = bit_helper_type::ALL_SET;
                }
            }
            
            process_space( ss, destF, B, pool );
        }

        this->analyze_free_indices( destF, destF + 2 * B, A );

        delete [] destF;
    }

    virtual ~FreeSpaceAnalyzerMT() {}

protected:

    template < class PoolType >
    void process_space( space_type * source, block_type * destF, const unsigned int BLOCK_COLUMNS, PoolType & pool ) {

        const unsigned int TC = pool.pool_size() + 1; // + 1 for master thread

        // batching by columns allows this algorithm to elimnate the need for mutex locks
        const unsigned int CPB = BLOCK_COLUMNS / TC + ((BLOCK_COLUMNS % TC > 0) ? 1 : 0);

        // cols will be decremented to CPB
        // unsigned int cols = BLOCK_COLUMNS;
        unsigned int cols = 0;
        while( cols + CPB < BLOCK_COLUMNS  ) {
            
            task_type t( source, destF, cols, cols + CPB );

            pool.post( t );
            cols += CPB;
        }

        if( cols < BLOCK_COLUMNS ) {
            // last block being run on master thread
            task_type t(source, destF,  cols, BLOCK_COLUMNS );
            t();
        }

        pool.sync();
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FREE_SPACE_MT_POPULATION_SPACE_ROW_HPP_

