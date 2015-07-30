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
#include <algorithm>

#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "clotho/cuda/curand_uniform_wrapper.hpp"
#include "clotho/cuda/curand_poisson_wrapper.hpp"

#include "crossover_matrix.hpp"
#include "population_recombiner.hpp"
#include "scatter_mutations.hpp"
#include "population_metadata.hpp"

#include <iostream>
#include <iomanip>

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"

#include <deque>

typedef clotho::utility::timer timer_type;

inline std::ostream & operator<<( std::ostream & out, const dim3 & d ) {
    out << "< " << d.x << ", " << d.y << ", " << d.z << " >";
    return out;
}

template < class Iter >
std::string to_string( Iter first, Iter last, std::string delim ) {
    if( first == last ) return "";

    std::ostringstream oss;
    oss << *first;
    ++first;
    while( first != last ) {
        if( !delim.empty() )
            oss << delim;

        oss << *first;
        ++first;
    }

    return oss.str();
}

template < class Iter >
std::string to_hex_string( Iter first, Iter last, std::string delim, unsigned int width = 8 ) {
    if( first == last ) return "";

    std::ostringstream oss;

    if( width ) {
        oss << std::hex << std::setfill('0') << std::setw(width) << *first;
        ++first;
        while( first != last ) {
            if( !delim.empty() )
                oss << delim;

            oss << std::hex << std::setfill('0') << std::setw(width) << *first;
            ++first;       
        }
    } else {
        oss << std::hex << *first;
        ++first;

        while( first != last ) {
            if( !delim.empty() )
                oss << delim;

            oss << std::hex << *first;
            ++first;
        }
    }

    return oss.str();
}

simulate_engine::simulate_engine( unsigned long seed, double mu, double rho, unsigned int founder_size ) :
    /*mutate_rate( mu )
    , recomb_rate( rho )
    ,*/
    m_dPop0()
    , m_dPop1()
    , m_dAlleleMasks()
    , m_dParentPop( &m_dPop0 )
    , m_dOffspringPop( &m_dPop1 )
    , m_dAlleles()
    , m_dOrderedAlleles()
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
    init( founder_size );
}

//simulate_engine::simulate_engine( const params & p ) {
//
//}

void simulate_engine::init( unsigned int founder_size ) {

    std::cerr << "Using Crossover Matrix Generator: " << CROSSOVER_MATRIX << std::endl;

    if( curandCreateGenerator( &m_dGen, CURAND_RNG_PSEUDO_MTGP32 ) != CURAND_STATUS_SUCCESS ) {
        assert(false);
    }

    typedef boost::random::uniform_int_distribution< unsigned long long > uniform_int_dist_type;
    uniform_int_dist_type uni;
    unsigned long long s = uni( m_hGen );

    std::cerr << "CURAND seed: " << s << std::endl;
    if( curandSetPseudoRandomGeneratorSeed( m_dGen, s ) != CURAND_STATUS_SUCCESS ) {
        assert(false);
    }

//    unsigned int a = 40000;
    unsigned int a = 1024;
    resizeAlleles(a);

    // initialize the alleles with random values
    clotho::cuda::fill_uniform< real_type > duni( m_dGen );
//    curand_gateway( m_dAlleles, m_dAlleles.size(), duni );
    duni( m_dAlleles, m_dAlleles.size() );

    thrust::fill(m_dFree.begin(), m_dFree.end(), -1);

    resizePopulation( m_dParentPop, 2 * founder_size );
}

void simulate_engine::simulate( unsigned int pop_size ) {
    swapPopulations();

    unsigned int seq_count = 2 * pop_size;  // diploid sequences

    unsigned int parent_alleles = m_dAlleles.size();
    unsigned int parent_cols = m_dFree.size();  // before resize of offspring

    unsigned int parent_rows = m_dParentPop->size() / parent_cols;

    assert( parent_rows * parent_cols == m_dParentPop->size() );

//    std::cerr << "Parent Alleles: " << parent_alleles << std::endl;
//    std::cerr << "Parent Dimensions: [" << parent_rows << ", " << parent_cols << "]" << std::endl;
    
//    timer_type t0;

    // pre-generate mutation event distribution for offspring population
    // done to allow offspring population to be resized appropriately
    unsigned int nMut = fill_event_list(m_dMutEvent, seq_count, m_mu);

    // resizeAlleles( nMut );
    resizePopulation( m_dOffspringPop, seq_count );

    // pre-generate recombination event distribution for offspring population
#if CROSSOVER_MATRIX == 2
    unsigned int nRec = fill_event_list(m_dRecEvent, seq_count, m_rho);
#else
    unsigned int nRec = fill_event_list( m_dRecEvent, seq_count * 1024, (m_rho / 1024.0));
#endif
//    std::cerr << "nRec: " << nRec << std::endl;

    // 1 (parent/offspring sequence) * (offspring sequences)
    // + mutations added to population
    // + recombination events
    fill_random_pool( seq_count + nMut + nRec );

//    t0.stop();
//    std::cerr << "Random Number Generation Lapse: " << t0 << std::endl;

//    t0.start();
    real_type * rand_pool = thrust::raw_pointer_cast( m_dRandBuffer.data() );
//
    crossover_method4( rand_pool, seq_count, parent_alleles );
    rand_pool += nRec;

    recombine_method3( rand_pool,  seq_count, parent_rows, parent_cols );
    rand_pool += seq_count;

    mutate_method1( rand_pool, seq_count, nMut );
    
//    update_population_metadata<<< blocks, threads >>>( offspring, seq_count, m_dFree.size(), free_list, lost_list, fixed_list );
    cudaDeviceSynchronize();
//    t0.stop();

//    std::cerr << "Generate Offspring Population Lapse: " << t0 << std::endl;
}

unsigned int simulate_engine::fill_event_list( event_vector & ev, unsigned int count, double rate ) {
    clotho::cuda::fill_poisson< unsigned int, double > event_gen( m_dGen, rate );
    curand_gateway( ev, count, event_gen );
    thrust::exclusive_scan( ev.begin(), ev.end(), ev.begin() );

    return ev.back();
}

void simulate_engine::fill_random_pool( size_t pool_size ) {
    clotho::cuda::fill_uniform< real_type > uni( m_dGen );
    curand_gateway( m_dRandBuffer, pool_size, uni );
}

void simulate_engine::crossover_method1( real_type * rand_pool, unsigned int seq_count, unsigned int nMut ) {

    // build recombination masks in child population space
    const unsigned int seq_per_block = 100;
    unsigned int bcount = seq_count / seq_per_block;
    if( seq_count % seq_per_block ) { ++bcount; }

    dim3 threads(32,32,1), blocks(m_dAlleles.size() / ALLELES_PER_STRIDE, bcount,1), sizes( m_dAlleles.size(), seq_count, 1 );

    real_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );
    event_count_type * rec_events = thrust::raw_pointer_cast( m_dRecEvent.data() );
    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );
    block_type * parents = thrust::raw_pointer_cast( m_dParentPop->data() );

    block_type * free_list = thrust::raw_pointer_cast( m_dFree.data() );

    // may need to specify a CUDA event
    generate_crossover_matrix<<< blocks, threads >>>( rand_pool, alleles, rec_events, offspring, sizes );
}

// inefficient. suspect too many serialized function calls
void simulate_engine::crossover_method2( unsigned int seq_count, unsigned int parent_alleles ) {

    // reset offspring population
    thrust::fill( m_dOffspringPop->begin(), m_dOffspringPop->end(), 0 );

    std::cerr << "Offspring size: " << m_dOffspringPop->size() << std::endl;

//    const unsigned int STREAM_COUNT = 8;
//    cudaStream_t streams[ STREAM_COUNT ];
//
//    std::deque< cudaStream_t * > incomplete( STREAM_COUNT, NULL );
//    
//    for( unsigned int i = 0; i < STREAM_COUNT; ++i ) {
//        cudaStreamCreate( &streams[i] );
//        incomplete[i] =  &streams[i];
//    }

    unsigned int bcount = m_dFree.size() / BLOCK_PER_STRIDE;
    if( m_dFree.size() % BLOCK_PER_STRIDE ) ++bcount;

    dim3 threads(BLOCK_PER_STRIDE,32,1), blocks( bcount, 1,1);

    boost::random::poisson_distribution< unsigned int, double > event_gen( m_rho );
    boost::random::uniform_real_distribution< double > rec_point_gen;

    unsigned int * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );
//    unsigned int * free_list = thrust::raw_pointer_cast( m_dFree.data() );
    double * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );

//    unsigned int sid = STREAM_COUNT;
    while(seq_count--) {
//        sid = (++sid % STREAM_COUNT);

        unsigned int n = event_gen( m_hGen );
        while( n-- ) {
            double r = rec_point_gen( m_hGen );
//            crossover<<< blocks, threads, streams[sid] >>>( offspring, alleles, free_list, parent_alleles, r );
//            crossover<<< blocks, threads, 0, streams[sid] >>>( offspring, alleles, parent_alleles, r );
            crossover<<< blocks, threads >>>( offspring, alleles, parent_alleles, r );
        }
        offspring += bcount;
    }

//    while( !incomplete.empty() ) {
//        cudaStream_t * t = incomplete.front();
//        incomplete.pop_front();
//
//        if( t == NULL ) continue;
//
//        if( incomplete.empty() || cudaStreamQuery( *t ) == cudaSuccess ) {
//            cudaStreamSynchronize( *t );
//
//            cudaStreamDestroy( *t );
//        } else {
//            incomplete.push_back( t );
//        }
//    }
    cudaDeviceSynchronize();
}

// inefficient. suspect too many serialized function calls
void simulate_engine::crossover_method3( unsigned int seq_count, unsigned int parent_alleles ) {

    allele_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );

    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );

    boost::random::poisson_distribution< unsigned int, double > event_gen( m_rho );
    boost::random::uniform_real_distribution< double > rec_point_gen;

    dim3 blocks(1,1,1), threads(32,32,1);
    while( parent_alleles ) {
        unsigned int N = parent_alleles;
        if( N > 1024 ) N = 1024;

        init_alleles<<< blocks, threads >>>( alleles, N );

        unsigned int i = 0;
        unsigned int * off = offspring;
        while( i < seq_count ) {
            unsigned int n = event_gen( m_hGen );
            init_sequence<<< blocks, threads >>>( off );

            while( n-- ) {
                double r = rec_point_gen( m_hGen );
                crossover2<<< blocks, threads >>>( r);
            }

            finalize_sequence<<< blocks, threads >>>( off );
            cudaDeviceSynchronize();
            off += m_dFree.size();
            ++i;
        }

        alleles += N;
        parent_alleles -= N;
        offspring += (N / 32);
    }
}

void simulate_engine::crossover_method4( real_type * rand_pool, unsigned int seq_count, unsigned int allele_count ) {
    const unsigned int STREAM_COUNT = 2;
    cudaStream_t streams[ STREAM_COUNT ];

    std::deque< cudaStream_t * > incomplete( STREAM_COUNT, NULL );
    
    for( unsigned int i = 0; i < STREAM_COUNT; ++i ) {
        cudaStreamCreate( &streams[i] );
        incomplete[i] =  &streams[i];
    }

    allele_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );
    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );
    event_count_type * event_list = thrust::raw_pointer_cast( m_dRecEvent.data() );

    unsigned int block_cols = (allele_count / ALLELES_PER_STRIDE );
    unsigned int max_block_rows = (comp_cap_type::MAX_BLOCKS_X / block_cols);
    dim3 threads( 32, 32, 1), sizes( m_dFree.size(), 0, 1);

    unsigned int sid = STREAM_COUNT;
    while( seq_count ) {
        sid = (++sid % STREAM_COUNT);

        dim3 blocks( block_cols, ((seq_count > max_block_rows) ? max_block_rows : seq_count) , 1 );

#if CROSSOVER_MATRIX == 2
        generate_crossover_matrix2a<<< blocks, threads, 0, streams[sid] >>>(rand_pool, alleles, event_list, offspring, sizes);
#else
        generate_crossover_matrix3<<< blocks, threads, 0, streams[sid] >>>(rand_pool, alleles, event_list, offspring, sizes);
#endif

        offspring += m_dFree.size() * blocks.y;
        seq_count -= blocks.y;
        sizes.y += blocks.y;
    }

    while( !incomplete.empty() ) {
        cudaStream_t * t = incomplete.front();
        incomplete.pop_front();

        if( t == NULL ) continue;

        if( incomplete.empty() || cudaStreamQuery( *t ) == cudaSuccess ) {
            cudaStreamSynchronize( *t );

            cudaStreamDestroy( *t );
        } else {
            incomplete.push_back( t );
        }
    }
}

void simulate_engine::recombine_method2( unsigned int seq_count, unsigned int parent_rows, unsigned int parent_cols ) {

    unsigned int offspring_cols = m_dFree.size();

    unsigned int * parents = thrust::raw_pointer_cast( m_dParentPop->data() );
    unsigned int * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );

    unsigned int bcount = parent_cols / 1024;
    if( parent_cols % 1024 ) {
        ++bcount;
    }

    dim3 blocks( bcount, 1, 1 ), threads(parent_cols, 1, 1);
    if( parent_cols > 1024 ) {
        threads.x = 1024;
    }

    boost::random::uniform_int_distribution< unsigned int > parent_gen( 0, (parent_rows / 2) - 1 );
    while( seq_count-- ) {
        unsigned int *p0 = parents + parent_gen( m_hGen ) * parent_cols;
        unsigned int *p1 = p0 + parent_cols;
        recombine<<< blocks, threads >>>( p0, p1, offspring, parent_cols );
        offspring += offspring_cols;
    }
}

void simulate_engine::recombine_method3( real_type * rand_pool, unsigned int seq_count, unsigned int p_row, unsigned int p_col ) {

    if( p_row == 0 ) {
        thrust::fill( m_dOffspringPop->begin(), m_dOffspringPop->end(), 0 );
        return;
    }

    unsigned int offspring_cols = m_dFree.size();

    unsigned int * parents = thrust::raw_pointer_cast( m_dParentPop->data() );
    unsigned int * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );

    // bad seed value: 4176498736
    dim3 blocks(1,1,1), threads(offspring_cols,1,1);

    if( offspring_cols > 1024) {
        threads.x = 1024;
        blocks.x = (offspring_cols / 1024) + 1;
    }

    unsigned int seq_offset = 0;
    while( seq_count ) {
        if( seq_count > comp_cap_type::MAX_BLOCKS_Y ) {
            blocks.y = comp_cap_type::MAX_BLOCKS_Y;
        } else {
            blocks.y = seq_count;
        }

        recombiner<<< blocks, threads >>>(rand_pool, parents, p_row, p_col, offspring, offspring_cols, seq_offset);

        seq_count -= blocks.y;
        offspring += (offspring_cols * blocks.y);
        seq_offset += blocks.y;
    }

    cudaDeviceSynchronize();
}

void simulate_engine::mutate_method1( real_type * rand_pool, unsigned int seq_count, unsigned int nMut ) {

    real_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );
    event_count_type * mut_events = thrust::raw_pointer_cast( m_dMutEvent.data());
    block_type * free_list = thrust::raw_pointer_cast( m_dFree.data() );
    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );

    dim3 blocks(1, 1, 1 ), threads(32, 32, 1);

    scatter_mutations<<< blocks, threads >>>( rand_pool, alleles, free_list, offspring, mut_events, seq_count, m_dAlleles.size() );
}

void simulate_engine::mutate_method2( real_type * rand_pool, unsigned int seq_count, unsigned int nMut ) {

    allele_type * alleles = thrust::raw_pointer_cast( m_dAlleles.data() );
    block_type * offspring = thrust::raw_pointer_cast( m_dOffspringPop->data() );
    block_type * free_list = thrust::raw_pointer_cast( m_dFree.data() );

    boost::random::poisson_distribution< unsigned int, double > event_gen( m_mu );
    boost::random::uniform_real_distribution< double > mut_gen( 0. ,1. );

    while( seq_count-- ) {
        unsigned int n = event_gen( m_hGen );
        while( n-- ) {
            double mut = mut_gen( m_hGen );

//            scatter<<< >>>( alleles, free_list, offspring, m_dFree.size(), mut );
        }
        offspring += m_dFree.size();
    }
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

//    std::cerr << "Allele Size: " << s << std::endl;

    m_dAlleles.resize(s);
    m_dAlleleMasks.resize(s);
    m_dOrderedAlleles.resize(s);

    s /= ALLELES_PER_BLOCK;

    m_dFree.resize( s );
    m_dFixed.resize( s );
    m_dLost.resize( s );

//    std::cerr << "resized size: " << m_dAlleles.size() << std::endl;
}

void simulate_engine::resizePopulation( data_vector * pop, size_t s ) {
//    std::cerr << "Resizing offspring: " << s * m_dFree.size() << std::endl;
    pop->resize( s * m_dFree.size() );
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
        out << " <<" << std::hex << std::setfill('0') << std::setw(8) << *first;
        while( ++first != last ) {
            if( i-- == 0 ) {
                out << ">\n<" << std::hex << std::setfill('0') << std::setw(8) << *first;
                i = block_width;
            } else {
                out << ", " << std::hex << std::setfill('0') << std::setw(8) << *first;
            }
        }
        out << ">>\n";
    } else {
        out << " Empty\n";
    }
}

template < class Iter >
void add_population_state( boost::property_tree::ptree & state, Iter first, Iter last, unsigned int n ) {
    while( first != last ) {
        clotho::utility::add_value_array( state, to_hex_string( first, first + n, ":" ) );
        first += n;
    }
}

void simulate_engine::record_state(  boost::property_tree::ptree & state ) {
    boost::property_tree::ptree r, m, u, a;
    clotho::utility::add_value_array(r, m_dRecEvent.begin(), m_dRecEvent.end() );
    state.add_child( "recombination.events", r );

    clotho::utility::add_value_array( m, m_dMutEvent.begin(), m_dMutEvent.end() );
    state.add_child( "mutations.events", m );

    clotho::utility::add_value_array( u, m_dRandBuffer.begin(), m_dRandBuffer.end() );
    state.add_child( "random.uniform_pool", u );

    clotho::utility::add_value_array( a, m_dAlleles.begin(), m_dAlleles.end() );
    state.add_child( "alleles.positions", a );

    boost::property_tree::ptree par, off;

    add_population_state( par, m_dParentPop->begin(), m_dParentPop->end(), m_dFree.size() );
    state.add_child( "parent.sequences", par);

    add_population_state( off, m_dOffspringPop->begin(), m_dOffspringPop->end(), m_dFree.size() );
    state.add_child( "offspring.sequences", off );

    state.put( "offspring.metadata.free_list", to_hex_string( m_dFree.begin(), m_dFree.end(), ":"));
    state.put( "offspring.metadata.fixed_list", to_hex_string(m_dFixed.begin(), m_dFixed.end(), ":"));
    state.put( "offspring.metadata.lost_list", to_hex_string(m_dLost.begin(), m_dLost.end(), ":"));
}

std::ostream & operator<<( std::ostream & out, const simulate_engine & se ) {
    out << "Simulate Engine:\n";

    out << "Alleles:";
    dump_vector( out, se.m_dAlleles.begin(), se.m_dAlleles.end() );

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

    out << "Free List:";
    dump_vector( out, se.m_dFree.begin(), se.m_dFree.end() );

    return out;
}
