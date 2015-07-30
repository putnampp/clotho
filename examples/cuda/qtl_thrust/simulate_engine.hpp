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
#ifndef SIMULATE_ENGINE_HPP_
#define SIMULATE_ENGINE_HPP_

#include <cuda.h>
#include <curand.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/property_tree/ptree.hpp>

#include <thrust/device_vector.h>

#include "clotho/cuda/compute_capability.hpp"

/*
class simulate_engine;
class test_log;
class test_log2;

class iStateLogger {
public:
    virtual void log_state( boost::property_tree::ptree & , const simulate_engine & ) = 0;
    virtual ~iStateLogger() {}
};

class mutate_rate {
public:
    friend class test_log;
    virtual ~mutate_rate() {}

protected:
    mutate_rate( double r ) : m_mu_r( r ) {}

private:
    double m_mu_r;
};

class recomb_rate {
public:
    friend class test_log2;
    virtual ~recomb_rate() {}

protected:
    recomb_rate( double r ) : m_rho_r( r ) {}

private:
    double m_rho_r;
};*/

class simulate_engine /*: public mutate_rate, public recomb_rate*/ {
public:
    typedef unsigned int    block_type;
    typedef double          allele_type;
    typedef double          real_type;
    typedef unsigned int    event_count_type;

    typedef compute_capability< 3, 0 > comp_cap_type;

    typedef thrust::device_vector< block_type >     data_vector;
    typedef thrust::device_vector< real_type >      random_real_vector;
    typedef thrust::device_vector< allele_type >    allele_vector;
    typedef thrust::device_vector< event_count_type >   event_vector;
    
    static const unsigned int WARP_SIZE = 32;
    static const unsigned int BLOCK_PER_WARP = 1;
    static const unsigned int BLOCK_PER_STRIDE = WARP_SIZE * BLOCK_PER_WARP; // 32
    static const unsigned int ALLELES_PER_BLOCK = sizeof( block_type ) * 8;  // 32
    static const unsigned int ALLELES_PER_STRIDE = BLOCK_PER_STRIDE * ALLELES_PER_BLOCK;// 32 * 32 = 1024

    simulate_engine( unsigned long seed, double mu, double rho, unsigned int founder );

    //simulate_engine( const params & p );

    void simulate( unsigned int pop_size );

    friend std::ostream & operator<<( std::ostream & out, const simulate_engine & se );

    void record_state( boost::property_tree::ptree & out );

    virtual ~simulate_engine();

protected:

    void init( unsigned int f);

    unsigned int fill_event_list(event_vector & ev, unsigned int n, double rate);
    void fill_random_pool( size_t pool_size );

    void crossover_method1( real_type * rand_pool, unsigned int seq_count, unsigned int nMut );
    void crossover_method2( unsigned int seq_count, unsigned int parent_alleles );
    void crossover_method3( unsigned int seq_count, unsigned int parent_alleles );
    void crossover_method4( real_type * rand_pool, unsigned int seq_count, unsigned int parent_alleles );

    void recombine_method2( unsigned int seq_count, unsigned int p_row, unsigned int p_col );
    void recombine_method3( real_type * rand_pool, unsigned int seq_count, unsigned int p_row, unsigned int p_col );

    void mutate_method1( real_type * rand_pool, unsigned int seq_count, unsigned int nMut );
    void mutate_method2( real_type * rand_pool, unsigned int seq_count, unsigned int nMut );

    template < class T, class Func >
    void    curand_gateway( thrust::device_vector< T > & buf, size_t N, Func & f ) {
        buf.resize( N + 1 );

        thrust::fill( buf.begin(), buf.end(), 0);

//        typedef boost::random::uniform_int_distribution< unsigned long long > uniform_int_dist_type;
//        uniform_int_dist_type uni;
//        unsigned long long seed = uni( m_hGen );
//
//        std::cerr << "CURAND seed: " << seed << std::endl;
//        if( curandSetPseudoRandomGeneratorSeed( m_dGen, seed ) != CURAND_STATUS_SUCCESS ) {
//        }

        f(buf, N);
    }

    void    pruneSpace();
    void    swapPopulations();

    void    resizeAlleles( size_t s );
    void    resizePopulation( data_vector * v, size_t s );


    data_vector         m_dPop0, m_dPop1, m_dAlleleMasks;
    data_vector         * m_dParentPop, * m_dOffspringPop;

    allele_vector       m_dAlleles, m_dOrderedAlleles;

    data_vector         m_dFree, m_dLost, m_dFixed;

    event_vector        m_dRecEvent, m_dMutEvent;
    random_real_vector  m_dRandBuffer;

    boost::random::mt19937  m_hGen;
    curandGenerator_t       m_dGen;

    double  m_mu, m_rho;
    size_t m_nFree;
};

std::ostream & operator<<( std::ostream & out, const simulate_engine & se );

/*
class test_log : public iStateLogger {
public:
    test_log( iStateLogger * o = NULL ) : m_dec(o) {}

    void log_state( boost::property_tree::ptree & p, const simulate_engine & se ) {
        if( m_dec ) {
            m_dec->log_state( p, se );
        }
        p.put("mu", se.m_mu_r );
    }
protected:
    iStateLogger * m_dec;
};

class test_log2 : public iStateLogger {
public:
    test_log2( iStateLogger * o = NULL ) : m_dec(o) {}

    void log_state( boost::property_tree::ptree & p, const simulate_engine & se ) {
        if( m_dec ) {
            m_dec->log_state( p, se );
        }
        p.put( "rho", se.m_rho_r );
        p.put( "mu2", se.m_mu_r );
    }
protected:
    iStateLogger * m_dec;
};*/

#endif  // SIMULATE_ENGINE_HPP_
