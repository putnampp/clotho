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

template < class SequenceSpaceType >
__global__ void all_samples( SequenceSpaceType * pop, basic_data_space< unsigned int > * samps ) {
    unsigned int N = pop->seq_count;

    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int tpb = blockDim.x * blockDim.y;

    if( tid == 0 ) {
        _resize_space_impl( samps, N );
    }
    __syncthreads();

    unsigned int * d = samps->data;
    while( tid < N ) {
        d[tid] = tid;
        tid += tpb;
    }
}


class PairwiseDifference : public clotho::utility::iStateObject {
public:

    PairwiseDifference( boost::property_tree::ptree & config ) {
        assert( cudaMalloc( (void **) &m_dStats, sizeof( pairwise_diff_stats ) ) == cudaSuccess);
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop, basic_data_space< unsigned int > * subpop ) {
        initialize_pairwise_diff_stats<<< 1, 32 >>>( m_dStats );
        CHECK_LAST_KERNEL_EXEC

//        dim3 threads(32, 32, 1);
        dim3 threads( 32, 1, 1);
        pairwise_difference_kernel<<< pairwise_diff_stats::BINS, threads >>>( pop->sequences.get_device_space(), subpop, m_dStats );
        CHECK_LAST_KERNEL_EXEC

        finalize_pairwise_diff_stats<<< 1, 32 >>>( m_dStats );
        CHECK_LAST_KERNEL_EXEC
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop ) {
        basic_data_space< unsigned int > * samps;
        create_space( samps );

        all_samples<<< 1, 1024 >>>( pop->sequences.get_device_space(), samps );

        evaluate( pop, samps );   

        delete_space( samps );
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree f;
        get_device_object_state( f, m_dStats );

        state.add_child( "difference", f );
    }

    virtual ~PairwiseDifference() {
        cudaFree( m_dStats );
    }

protected:
    pairwise_diff_stats     * m_dStats;
};

#endif  // CUDA_PAIRWISE_DIFFERENCE_HPP_
