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

#include "clotho/cuda/data_spaces/data_space.hpp"

#include "clotho/cuda/analysis/sequence_hamming_weight_kernel.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/popcount.hpp"

#include <boost/random/mersenne_twister.hpp>

typedef unsigned int int_type;

template < unsigned char V >
struct test_kernel {
    test_kernel() {
        create_space( dWeights );
    }

    void operator()( device_sequence_space< int_type > * seqs, unsigned int N ) {
        resize_space( dWeights, N );

        algo_version< V > * v = NULL;
        sequence_hamming_weight_kernel<<< 10, 1024 >>>( seqs, dWeights, v );
    }

    void validate( int_type * exp_weights, unsigned int N ) {
        basic_data_space< int_type > obs;
        assert( cudaMemcpy( &obs, dWeights, sizeof( basic_data_space< int_type > ), cudaMemcpyDeviceToHost ) == cudaSuccess );

        assert( obs.size == N );

        int_type * obs_weights = new int_type[ N ];

        copy_heap_data(obs_weights, obs.data, N );

        unsigned int i = 0;

        while( i < N ) {
            if( exp_weights[ i ] != obs_weights[ i ] ) {
                std::cerr << "ERROR -";
                std::cerr << "\nKernel version = " << (unsigned int) V;
                std::cerr << "\nMessage = Weights do not match";
                std::cerr << "\nIndex = " << i;
                std::cerr << "\nObserved Weight = " << obs_weights[ i ];
                std::cerr << "\nExpected Weight = " << exp_weights[ i ];
                std::cerr << std::endl;
    
                break;
            }
            ++i;
        }

        delete obs_weights;
    }

    virtual ~test_kernel() {
        delete_space( dWeights );
    }

    basic_data_space< int_type > * dWeights;
};

void compute_weights( device_sequence_space< int_type > & seqs, int_type * weights ) {
    unsigned int cols = seqs.seq_width;
    unsigned int N = seqs.size;

    unsigned int idx = 0, count = 0;
    while( idx < N ) {
        int_type b = seqs.sequences[ idx++ ];
        count += popcount( b );

        if( idx % cols == 0 ) {
            weights[ (idx / cols) - 1 ] = count;
            count = 0;
        }
    }
}

void upload_sequences( device_sequence_space< int_type > & seqs, device_sequence_space< int_type > *dSeqs ) {
    device_sequence_space< int_type > local;

    assert( cudaMemcpy( &local, dSeqs, sizeof( device_sequence_space< int_type > ), cudaMemcpyDeviceToHost) == cudaSuccess );

    unsigned int N = seqs.size;
    assert( local.size == N );

    int_type * dTmp;

    N *= sizeof( int_type );
    assert( cudaMalloc( (void **) &dTmp, N ) == cudaSuccess );

    assert( cudaMemcpy( dTmp, seqs.sequences, N, cudaMemcpyHostToDevice ) == cudaSuccess );

    copy_heap<<< 1, 1 >>>( local.sequences, dTmp, local.size );
    cudaDeviceSynchronize();

    cudaFree( dTmp );
}

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

int main( int argc, char ** argv ) {

    device_sequence_space< int_type >   seqs;
    seqs.sequences = NULL;

    build_sequence_space( &seqs );

    device_sequence_space< int_type >   * dSeqs;
    create_space( dSeqs );

    resize_space( dSeqs, seqs.seq_width, seqs.seq_count );
    upload_sequences( seqs, dSeqs );

    test_kernel< 1 > t1;
    test_kernel< 2 > t2;
    test_kernel< 3 > t3;
    test_kernel< 4 > t4;

    t1( dSeqs, seqs.seq_count );
    t2( dSeqs, seqs.seq_count );
    t3( dSeqs, seqs.seq_count );
    t4( dSeqs, seqs.seq_count );

    int_type    * exp_weights = new int_type[ seqs.seq_count ];
    compute_weights( seqs, exp_weights );

    t1.validate( exp_weights, seqs.seq_count );
    t2.validate( exp_weights, seqs.seq_count );
    t3.validate( exp_weights, seqs.seq_count );
    t4.validate( exp_weights, seqs.seq_count );

    unsigned int N = ((10 < seqs.seq_count ) ? 10 : seqs.seq_count );

    std::cout << "Sample of Expected Weights: < " << exp_weights[0];
    for( unsigned int i = 1; i < N; ++i ) {
        std::cout << ", " << exp_weights[ i ];
    }
    std::cout << ", ... >" << std::endl;

    if( seqs.sequences ) {
        delete seqs.sequences;
    }

    delete exp_weights;

    delete_space( dSeqs );

    return 0;
}
