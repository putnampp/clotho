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
#include "clotho/cuda/execution_configuration/thread_parameter.hpp"
#include "clotho/cuda/execution_configuration/block_parameter.hpp"

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

    PairwiseDifference( boost::property_tree::ptree & config ) :
        m_dStats(NULL)
        , m_threads( 32, 32, 1)
        , m_kern_version( 0 )
    {
        assert( cudaMalloc( (void **) &m_dStats, sizeof( pairwise_diff_stats ) ) == cudaSuccess);

        parse_configuration( config );
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop, basic_data_space< unsigned int > * subpop ) {
        initialize_pairwise_diff_stats<<< 1, 32 >>>( m_dStats );
        CHECK_LAST_KERNEL_EXEC

        if( m_kern_version ) {
            pairwise_difference_kernel2<<< pairwise_diff_stats::BINS, m_threads >>>( pop->sequences.get_device_space(), subpop, m_dStats );
        } else {
            pairwise_difference_kernel<<< pairwise_diff_stats::BINS, m_threads >>>( pop->sequences.get_device_space(), subpop, m_dStats );
        }
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
    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lc, kern, th;

        lc = config.get_child( "pairwise_difference", lc);

        kern = lc.get_child( "kernel", kern );

        m_kern_version = kern.get< unsigned int >( "version", m_kern_version );
        kern.put( "version", m_kern_version );

        th = kern.get_child( "threads", th );

        m_threads.x = th.get< unsigned int >( "x", m_threads.x );
        m_threads.y = th.get< unsigned int >( "y", m_threads.y );
        m_threads.z = th.get< unsigned int >( "z", m_threads.z );

        th.put( "x", m_threads.x);
        th.put( "y", m_threads.y);
        th.put( "z", m_threads.z);

        kern.put_child( "threads", th );

        lc.put_child( "kernel", kern );   

        config.put_child( "pairwise_difference", lc );
    }

    pairwise_diff_stats     * m_dStats;
    dim3 m_threads;
    unsigned int m_kern_version;
};

#endif  // CUDA_PAIRWISE_DIFFERENCE_HPP_
