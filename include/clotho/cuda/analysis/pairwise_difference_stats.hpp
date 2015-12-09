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
#ifndef CUDA_PAIRWISE_DIFFERENCE_STATS_HPP_
#define CUDA_PAIRWISE_DIFFERENCE_STATS_HPP_

#include "clotho/utility/state_object.hpp"

/// TODO need to figure out how to compute standard deviation
struct pairwise_diff_stats {
    static const unsigned int BINS = 32;
    
    unsigned long long block_bin[ BINS ];

    unsigned long long count, total;
    double mean;
};

namespace clotho {
namespace utility {

template <>
struct state_getter< pairwise_diff_stats > {
    void operator()( boost::property_tree::ptree & state, const pairwise_diff_stats & stats ) {
        state.put( "mean", stats.mean );
        state.put( "count", stats.count );
        state.put( "total", stats.total );

    }
};

}
}

#endif  // CUDA_PAIRWISE_DIFFERENCE_STATS_HPP_
