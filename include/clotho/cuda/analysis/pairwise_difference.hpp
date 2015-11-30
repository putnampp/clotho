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
#ifndef CUDA_PAIRWISE_DIFFERENCE_HPP_
#define CUDA_PAIRWISE_DIFFERENCE_HPP_

#include "clotho/cuda/device_state_object.hpp"

#include "clotho/cuda/analysis/pairwise_difference_kernel.hpp"

#include "clotho/cuda/helper_macros.hpp"

class PairwiseDifference : public clotho::utility::iStateObject {
public:

    PairwiseDifference( boost::property_tree::ptree & config ) {
        assert( cudaMalloc( (void **) &m_dStats, sizeof( pairwise_diff_stats ) ) == cudaSuccess);
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop ) {
        dim3 threads(32, 32, 1);
        pairwise_difference_kernel<<< pairwise_diff_stats::BINS, threads >>>( pop->sequences.get_device_space(), m_dStats );
        CHECK_LAST_KERNEL_EXEC
    }

    virtual ~PairwiseDifference() {
        cudaFree( m_dStats );
    }

protected:
    pairwise_diff_stats     * m_dStats;
};

#endif  // CUDA_PAIRWISE_DIFFERENCE_HPP_
