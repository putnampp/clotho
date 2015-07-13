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

#include <thrust/device_vector.h>

class simulate_engine {
public:
    typedef unsigned int    block_type;
    typedef double          allele_type;
    typedef double          real_type;
    typedef unsigned int    event_count_type;

    typedef thrust::device_vector< block_type >     data_vector;
    typedef thrust::device_vector< real_type >      random_real_vector;
    typedef thrust::device_vector< allele_type >    allele_vector;
    typedef thrust::device_vector< event_count_type >   event_vector;
    
    static const unsigned int WARP_SIZE = 32;
    static const unsigned int BLOCK_PER_WARP = 1;
    static const unsigned int BLOCK_PER_STRIDE = WARP_SIZE * BLOCK_PER_WARP; // 32
    static const unsigned int ALLELES_PER_BLOCK = sizeof( block_type ) * 8;  // 32
    static const unsigned int ALLELES_PER_STRIDE = BLOCK_PER_STRIDE * ALLELES_PER_BLOCK;// 32 * 32 = 1024

    simulate_engine( unsigned long seed, double mu, double rho );

    //simulate_engine( const params & p );

    void simulate( unsigned int pop_size );

    friend std::ostream & operator<<( std::ostream & out, const simulate_engine & se );

    virtual ~simulate_engine();

protected:

    template < class T, class Func >
    void    curand_gateway( thrust::device_vector< T > & buf, size_t N, Func & f ) {
        buf.resize( N + 1 );

        thrust::fill( buf.begin(), buf.end(), 0);

        typedef boost::random::uniform_int_distribution< unsigned long long > uniform_int_dist_type;
        uniform_int_dist_type uni;
        unsigned long long seed = uni( m_hGen );

        std::cerr << "CURAND seed: " << seed << std::endl;
        if( curandSetPseudoRandomGeneratorSeed( m_dGen, seed ) != CURAND_STATUS_SUCCESS ) {
        }

        f(buf, N);
    }

    void    pruneSpace();
    void    swapPopulations();

    void    resizeAlleles( size_t s );
    void    resizeOffspring( size_t s );


    data_vector         m_dPop0, m_dPop1;
    data_vector         * m_dParentPop, * m_dOffspringPop;

    allele_vector       m_dAlleles;

    data_vector         m_dFree, m_dLost, m_dFixed;

    event_vector        m_dRecEvent, m_dMutEvent;
    random_real_vector  m_dRandBuffer;

    boost::random::mt19937  m_hGen;
    curandGenerator_t       m_dGen;

    double  m_mu, m_rho;
    size_t m_nFree;
};

std::ostream & operator<<( std::ostream & out, const simulate_engine & se );

#endif  // SIMULATE_ENGINE_HPP_
