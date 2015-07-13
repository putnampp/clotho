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
#include "simulate_engine.hpp"

#include <algorithm>
#include <thrust/scan.h>

#include "curand_uniform_wrapper.hpp"
#include "curand_poisson_wrapper.hpp"

#include "crossover_matrix.hpp"

#include <iostream>

#include "clotho/utility/timer.hpp"

const unsigned int THREAD_COUNT = 1024;

//__global__ void updateReproductionPage( unsigned int ** repro_pages, unsigned int ** page_digest )

/**
 * Uncertain whether this version will still access main memory in strides
 * Syntactically simpler code
 */
//__global__ void computeParents( unsigned int * parents, float * offsets, unsigned int M, unsigned int N ) {
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//
//    unsigned int idx = bid * (blockDim.x * blockDim.y) + tid;
//
//    if( idx < N ) {
//        parents[ idx ] = offsets[idx] * (M - 1);
//    }
//}

//__global__ void computeParentsShared( unsigned int * parents, float * offsets, unsigned int M, unsigned int N ) {
//
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//
//    unsigned int idx = bid * (blockDim.x * blockDim.y) + tid;
//
//    __shared__ float loff[ THREAD_COUNT ];
//
//    loff[ tid ] = ((idx < N) ? offsets[idx] : 0f);
//    __syncthreads();
//    
//    if( idx < N ) {
//        parents[ idx ] = loff[tid] * (M - 1);
//    }
//}

//__global__ void reproduceParents( unsigned int * parent, unsigned int * child ) {
//
//    unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
//    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
//
//    unsigned int offset = (bid * (blockDim.y) + threadIdx.y) / 2;   // Does this perform integer division for offset correctly? Think so.
//
//    __shared__ unsigned int buffer[ THREAD_COUNT ];
//
//    // broadcast parent offset to all threads
//    buffer[tid] = parent[ offset ];
//    __syncthreads();
//
//    // broadcast sequence_offset to all threads in parent warps
//    offset = 2 * buffer[tid] + (threadIdx.y & 1);
//    buffer[tid] = sequence_lookup[ offset ];
//    __syncthreads();
//
//    offset = buffer[tid] + threadIdx.x;
//    buffer[tid] = sequences[ offset ];
//}

simulate_engine::simulate_engine( unsigned long seed, double mu, double rho ) :
    m_dPop0()
    , m_dPop1()
    , m_dParentPop( &m_dPop0 )
    , m_dOffspringPop( &m_dPop1 )
    , m_dAlleles()
    , m_dFree()
    , m_dLost()
    , m_dFixed()
    , m_dRecEvent()
    , m_dRandBuffer()
    , m_hGen( seed )
    , m_dGen(NULL)
    , m_mu(mu)
    , m_rho(rho)
    , m_nFree(0)
{
    if( curandCreateGenerator( &m_dGen, CURAND_RNG_PSEUDO_MTGP32 ) != CURAND_STATUS_SUCCESS ) {

    }
}

//simulate_engine::simulate_engine( const params & p ) {
//
//}

void simulate_engine::simulate( unsigned int pop_size ) {

    unsigned int seq_count = 2 * pop_size;  // diploid sequences

    typedef clotho::utility::timer timer_type;

    timer_type t0;

    fill_poisson< unsigned int, double > mut_events_gen( m_dGen, m_mu ), rec_events_gen( m_dGen, m_rho );

    // generate number of mutation and recombination events for offspring population
    curand_gateway( m_dMutEvent, seq_count, mut_events_gen );
    curand_gateway( m_dRecEvent, seq_count, rec_events_gen );

    thrust::exclusive_scan( m_dMutEvent.begin(), m_dMutEvent.end(), m_dMutEvent.begin() );
    thrust::exclusive_scan( m_dRecEvent.begin(), m_dRecEvent.end(), m_dRecEvent.begin() );

    unsigned int nMut = m_dMutEvent.back();
    unsigned int nRec = m_dRecEvent.back();

    t0.stop();

    std::cerr << "< " << nMut << ", " << nRec << " > ->  < " << m_dMutEvent.size() << ", " << m_dRecEvent.size() << "> -> " << t0 << std::endl;

    resizeAlleles( nMut );

    t0.start();
    // resize offspring population (if necessary)
    resizeOffspring( seq_count );

    // seq_count = 2 * pop_size = # of parents
    // nMut = # of mutations to be generated
    // nRec = # of recombination events to generate
    fill_uniform< real_type > uni( m_dGen );

    // child sequences + nMut + nRec
    unsigned int nVariables = seq_count + nMut + nRec;
    curand_gateway( m_dRandBuffer, nVariables, uni );

    t0.stop();

    std::cerr << "Uniform variables: " << nVariables << " -> " << t0 << std::endl;

    // build recombination masks in child population space
    dim3 threads(32,32,1), blocks(m_dAlleles.size() / ALLELES_PER_STRIDE,1,1), sizes( m_dAlleles.size(), seq_count, nVariables );

    real_type * unirands = thrust::raw_pointer_cast( m_dRandBuffer.data() );
    real_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );
    event_count_type * rec_events = thrust::raw_pointer_cast( m_dRecEvent.data() );
    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );
    block_type * parents = thrust::raw_pointer_cast( m_dParentPop->data() );

    // may need to specify a CUDA event
    generate_crossover_matrix<<< blocks, threads >>>( unirands, rec_events, alleles, offspring, sizes );

    unirands += nRec;   // move pointer past the recombination events

    //recombine<<< blocks, threads >>>( unirand, parents, offspring, sizes );

    // scatter new mutations in child population sequences in device memory
    //scatterMutationInOffspring<<<grid_size, block_size>>>( m_dOffspringPop, pop_size, nMut, m_dFreeAlleles );

//    swapPopulations();
//    pruneSpace();
}

void simulate_engine::swapPopulations() {
    std::swap( m_dParentPop, m_dOffspringPop );
}

void simulate_engine::pruneSpace() {

}

void simulate_engine::resizeAlleles( size_t s ) {
    // pad s to be a multiple of ALLELES_PER_STRIDE

    if( s < m_nFree ) { return; }

    s += m_dAlleles.size() - m_nFree;

    unsigned int tail = s % ALLELES_PER_STRIDE;
    if( tail ) { s += (ALLELES_PER_STRIDE - tail); }

    m_dAlleles.resize(s);

    s /= ALLELES_PER_BLOCK;

    m_dFree.resize( s );
    m_dFixed.resize( s );
    m_dLost.resize( s );
}

void simulate_engine::resizeOffspring( size_t s ) {
    m_dOffspringPop->resize( s * m_dFree.size() );
}

simulate_engine::~simulate_engine() {
}

template < class Iter >
void dump_vector( std::ostream & out, Iter first, Iter last ) {
    if( first != last ) {
        out << " <" << *first;
        while( ++first != last ) {
            out << ", " << *first;
        }
        out << ">\n";
    } else {
        out << " Empty\n";
    }
}

template < class Iter >
void dump_vector( std::ostream & out, Iter first, Iter last, unsigned int block_width ) {
    if( first != last ) {
        unsigned int i = --block_width;
        out << " <<" << *first;
        while( ++first != last ) {
            if( i-- == 0 ) {
                out << ">\n<" << *first;
                i = block_width;
            } else {
                out << ", " << *first;
            }
        }
        out << ">>\n";
    } else {
        out << " Empty\n";
    }
}

std::ostream & operator<<( std::ostream & out, const simulate_engine & se ) {
    out << "Simulate Engine:\n";

    out << "Recombination Events:";
    dump_vector( out, se.m_dRecEvent.begin(), se.m_dRecEvent.end(), 100 );

    out << "Mutation Events:";
    dump_vector( out, se.m_dMutEvent.begin(), se.m_dMutEvent.end(), 100 );

    out << "Uniform Numbers:";
    dump_vector( out, se.m_dRandBuffer.begin(), se.m_dRandBuffer.end() );

    out << "Child Population:";
    dump_vector( out, se.m_dOffspringPop->begin(), se.m_dOffspringPop->end(), se.m_dFree.size() );

    out << "Parent Population:";
    dump_vector( out, se.m_dParentPop->begin(), se.m_dParentPop->end(), se.m_dFree.size() );
    return out;
}
