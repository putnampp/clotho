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
#ifndef CLOTHO_CROSSOVER_TASK_LIST_HPP_
#define CLOTHO_CROSSOVER_TASK_LIST_HPP_

#include "clotho/data_spaces/crossover/tail_crossover_task.hpp"
#include "clotho/data_spaces/crossover/segment_crossover_task.hpp"
#include "clotho/data_spaces/crossover/copy_crossover_task.hpp"

#include <boost/random/bernoulli_distribution.hpp>

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType, class RNG >
static task_sequence make_crossover_task_sequence( const Classifier & events, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, double strand_bias, RNG & rand ) {

    boost::random::bernoulli_distribution< double > dist( strand_bias );

    if( dist(rand) ) {
        // consider s1 as top strand, s0 as bottom strand
        std::swap( s0, s1 );
        std::swap( s0_len, s1_len );
    }

    TaskGroup tg;    

    if( events.count() == 0 ) {
        // there are no crossover events
        // therefore, offspring strand will be a copy of the top strand
        BlockType * res = new BlockType[];
        tg.add( copy_crossover_task( events, s1, res, s1_len ) );
    } else if( s0_len < s1_len ) {
        // top strand is shorter than bottom strand
        BlockType res = new BlockType[];
        tg.add( segment_crossover_task( events, s0, s1, 0, s0_len ) );
        tg.add( tail_crossover_task( events, s1, res, s0_len, s1_len ) );
    } else if( s1_len < s0_len ) {
        // bottom strand is shorter than top strand
        BlockType res = new BlockType[ ];
        tg.add( segment_crossover_task( events, s0, s1, 0, s1_len ) );
        tg.add( tail_crossover_task( events, s0, res, s1_len, s0_len ) );
    } else {
        // both strands are equal in length
        BlockType res = new BlockType[ ];
        tg.add( segment_crossover_task( events, s0, s1, res, 0, s0_len ) );
    }

    return tg;
}


}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_TASK_LIST_HPP_
