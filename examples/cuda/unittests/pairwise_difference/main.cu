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
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

#include "clotho/cuda/data_spaces/data_space.hpp"

#include "clotho/cuda/analysis/pairwise_difference_kernel.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/popcount.hpp"
#include "clotho/utility/algorithm_version.hpp"

#include "clotho/genetics/pairwise_statistic.hpp"

#include "clotho/powerset/variable_subset.hpp"

#include <boost/random/mersenne_twister.hpp>

typedef unsigned int int_type;

typedef clotho::powersets::variable_subset< int, int_type > subset_type;
typedef typename subset_type::powerset_type powerset_type;

template < class IntType >
__global__ void create_subpopulation( basic_data_space< IntType > * data, unsigned int N ) {
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;

    unsigned int tpb = blockDim.x * blockDim.y;

    if( tid == 0 ) {
        _resize_space_impl( data, N );
    }
    __syncthreads();

    IntType * d = data->data;

    while( tid < N ) {
        d[tid] = tid;

        tid += tpb;        
    }
}

struct test_kernel {
    test_kernel() {
        create_space( dSubPop );
        assert( cudaMalloc( (void **)&dStats, sizeof(pairwise_diff_stats) ) == cudaSuccess);
    }

    void operator()( device_sequence_space< int_type > * seqs, unsigned int N ) {
        create_subpopulation<<< 1, 1024 >>>( dSubPop, N );

        initialize_pairwise_diff_stats<<< 1, 32 >>>( dStats );

        pairwise_difference_kernel<<< 10, 32 >>>( seqs, dSubPop, dStats);

        finalize_pairwise_diff_stats<<< 1, 32 >>>( dStats );
    }

    void validate( pairwise_diff_stats & exp_stats ) {
        pairwise_diff_stats obs_stats;

        assert( cudaMemcpy( &obs_stats, dStats, sizeof(pairwise_diff_stats), cudaMemcpyDeviceToHost ) == cudaSuccess );

        if( abs(exp_stats.mean - obs_stats.mean) > 0.00001 ) {
            std::cerr << "ERROR -";
            std::cerr << "\nMessage = Average pairwise difference does not match";
            std::cerr << "\nObserved Mean = " << obs_stats.mean;
            std::cerr << "\nObserved Count = " << obs_stats.count;
            std::cerr << "\nObserved Total = " << obs_stats.total;
            std::cerr << "\nExpected Mean = " << exp_stats.mean;
            std::cerr << "\nExpected Count = " << exp_stats.count;
            std::cerr << std::endl;
        } else {
            std::cout << "Expected Mean = " << exp_stats.mean;
            std::cout << "\nObserved Mean = " << obs_stats.mean;
            std::cout << std::endl;
        }
    }

    virtual ~test_kernel() {
        delete_space( dSubPop );
        cudaFree(dStats);
    }

    basic_data_space< int_type > *  dSubPop;
    pairwise_diff_stats         *  dStats;
};

void upload_sequences( device_sequence_space< int_type > & seqs, device_sequence_space< int_type > *dSeqs ) {
    device_sequence_space< int_type > local;

    assert( cudaMemcpy( &local, dSeqs, sizeof( device_sequence_space< int_type > ), cudaMemcpyDeviceToHost) == cudaSuccess );

    unsigned int N = seqs.size;
    assert( local.size == N );

    int_type * dTmp;

    N *= sizeof( int_type );
    assert( cudaMalloc( (void **) &dTmp, N ) == cudaSuccess );

    assert( cudaMemcpy( dTmp, seqs.sequences, N, cudaMemcpyHostToDevice ) == cudaSuccess );

    copy_heap<<< 1, 1024 >>>( local.sequences, dTmp, local.size );
    cudaDeviceSynchronize();

    cudaFree( dTmp );
}

/*
 * NOTE: device_sequence_space and basic_data_space objects allocate memory in the device's HEAP SPACE
 * The size of the device's HEAP SPACE is configurable.  It's default value is 8MB.
 * If expanding max_dim > 200, then consider increasing the size of the HEAP SPACE
 */
void build_sequence_space( device_sequence_space< int_type > * seqs, unsigned int max_dim = 200) {
    unsigned long long seed = clotho::utility::clock_type::now().time_since_epoch().count();
    boost::random::mt19937  rand_eng( seed );

    seqs->seq_width = (rand_eng() % max_dim) + 1;
    seqs->seq_count = (rand_eng() % max_dim) + 1;

    seqs->size = seqs->seq_width * seqs->seq_count;
    seqs->capacity = seqs->size;

    std::cout << "Constructing sequence space: < " << seqs->seq_count << ", " << seqs->seq_width << " >" << std::endl;

    if( seqs->sequences ) {
        delete seqs->sequences;
    }

    seqs->sequences = new int_type[ seqs->size ];

    unsigned int i = 0;
    while ( i < seqs->size ) {
        seqs->sequences[ i++ ] = rand_eng();    // fill sequence with random bits
    }
}

void compute_expectation( device_sequence_space< int_type > &  seqs, pairwise_diff_stats & exp_stats ) {
    powerset_type pset;

    typedef typename powerset_type::dc_type dc_type;
    typedef typename powerset_type::subset_ptr subset_ptr;

    typedef std::vector< subset_ptr >   population_type;

    population_type pop;

    unsigned int offset = 0;
    for( unsigned int i = 0; i < seqs.seq_count; ++i ) {
        dc_type cont;
        for( unsigned int j = 0; j < seqs.seq_width; ++j ) {
            cont.append( seqs.sequences[offset++] );
        }

        pop.push_back( pset.create_subset( cont ) );
    }

    unsigned int N = seqs.seq_count;
    exp_stats.count = (N * (N - 1) ) / 2;
    typename pairwise_statistic< subset_type >::accum_type global_diff, global_int, global_un;

    typename population_type::iterator kit = pop.begin();
    for( ; kit != pop.end(); ++kit ) {
        pairwise_statistic< subset_type > pstat( **kit, global_diff, global_int, global_un, 1.0 );
        for( typename population_type::iterator  kit2 = kit + 1; kit2 != pop.end(); ++kit2 ) {
            pstat.update( **kit2, 1.0 );
        }
    }

    exp_stats.mean = accum::weighted_mean(global_diff);
}

int main( int argc, char ** argv ) {

    device_sequence_space< int_type >   seqs;
    seqs.sequences = NULL;

    build_sequence_space( &seqs );

    device_sequence_space< int_type >   * dSeqs;
    create_space( dSeqs );

    resize_space( dSeqs, seqs.seq_width, seqs.seq_count );
    upload_sequences( seqs, dSeqs );

    test_kernel t1;

    t1( dSeqs, seqs.seq_count );

    pairwise_diff_stats exp_stats;
    compute_expectation( seqs, exp_stats );

    t1.validate( exp_stats );

    if( seqs.sequences ) {
        delete seqs.sequences;
    }

    delete_space( dSeqs );

    return 0;
}
