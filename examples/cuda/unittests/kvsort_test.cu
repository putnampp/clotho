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
#include <cassert>

#include <cuda.h>
#include <curand.h>
#include <thrust/device_vector.h>

#include <boost/lexical_cast.hpp>

#include "clotho/cuda/compute_capability.hpp"
#include "clotho/cuda/order_warp/order_warp.hpp"
#include "clotho/cuda/curand_uniform_wrapper.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/utility/log_helper.hpp"

#include <map>

typedef double key_type;
typedef unsigned int value_type;
typedef compute_capability< 3, 0 > comp_cap_type;

typedef clotho::cuda::order_warp< comp_cap_type > sort_type;

template < class K, class V >
void periodic_sort( thrust::device_vector< K > & kvec
                    , thrust::device_vector< V > & vvec
                    , thrust::device_vector< K > & sorted_keys
                    , thrust::device_vector< V > & sorted_values ) {

    assert( kvec.size() == vvec.size() );

    unsigned int N = kvec.size();
    sort_type sorter;

    sorted_keys.resize( N );
    sorted_values.resize( N );

    K * keys = thrust::raw_pointer_cast( kvec.data() );
    V * vals = thrust::raw_pointer_cast( vvec.data() );

    K * skeys = thrust::raw_pointer_cast( sorted_keys.data() );
    V * svals = thrust::raw_pointer_cast( sorted_values.data() );
    sorter.sort( keys, vals, skeys, svals, N );
}

template < class K, class V >
void dump( thrust::device_vector< K > & keys, thrust::device_vector< V > & values, boost::property_tree::ptree & log ) {
    assert( keys.size() == values.size() );

    typedef typename thrust::device_vector< K >::iterator key_iterator;
    typedef typename thrust::device_vector< V >::iterator value_iterator;

    key_iterator kit = keys.begin();
    value_iterator vit = values.begin();
    while( kit != keys.end() ) {
        clotho::utility::add_value_array( log, std::make_pair( *kit, *vit ) );
        ++kit;
        ++vit;
    }
}

template < class K, class V >
bool validate_sorted(thrust::device_vector< K > & kvec
                    , thrust::device_vector< V > & vvec
                    , thrust::device_vector< K > & sorted_keys
                    , thrust::device_vector< V > & sorted_values
                    , boost::property_tree::ptree & err ) {
    bool success = (kvec.size() == vvec.size() && sorted_keys.size() == sorted_values.size() && kvec.size() == sorted_keys.size());

    typedef std::map< K, V > map_type;
    typedef typename map_type::iterator pair_iterator;

    typedef typename thrust::device_vector< K >::iterator key_iterator;
    typedef typename thrust::device_vector< V >::iterator value_iterator;

    key_iterator k_it = kvec.begin();
    key_iterator k_sit = sorted_keys.begin();
    value_iterator v_it = vvec.begin();
    value_iterator v_sit = sorted_values.begin();

    unsigned int pair_idx = 0;

    while( success ) {
        map_type    period;

        unsigned int i = comp_cap_type::WARP_SIZE;
        while( i-- && k_it != kvec.end() ) {
            K k = *k_it;
            V v = *v_it;

            period.insert( std::make_pair( k, v ) );
            ++k_it;
            ++v_it;
        }

        pair_iterator p_it = period.begin();
        while(success && p_it != period.end() ) {
            success = success && (p_it->first == *k_sit) && (p_it->second == *v_sit);

            if( !success ) {
                err.put("error.pair_index", pair_idx );
                err.put("error.expected.key", p_it->first );
                err.put("error.expected.value", p_it->second );
                err.put("error.observed.key", *k_sit);
                err.put("error.observed.value", *v_sit);
            }
            ++p_it;
            ++k_sit;
            ++v_sit;
            ++pair_idx;
        }

        if( k_it == kvec.end() || k_sit == sorted_keys.end() || v_it == vvec.end() || v_sit == sorted_values.end() ) {
            success = success && (k_it == kvec.end() && k_sit == sorted_keys.end() && v_it == vvec.end() && v_sit == sorted_values.end());
            break;
        }
    }

    return success;
}


int main( int argc, char ** argv ) {

    if( argc != 3 ) {
        std::cerr << "<program> <seed> <N>" << std::endl;
        return 1;
    }

    unsigned long long seed = boost::lexical_cast< unsigned long long >( argv[1] );
    unsigned int N = boost::lexical_cast< unsigned int >( argv[2] );

    thrust::device_vector< double > key_vec, sorted_keys;
    thrust::device_vector< unsigned int > val_vec, sorted_vals;

    curandGenerator_t dGen;
    if( curandCreateGenerator( &dGen, CURAND_RNG_PSEUDO_MTGP32 ) != CURAND_STATUS_SUCCESS ) {
        assert(false);
    }

    if( curandSetPseudoRandomGeneratorSeed( dGen, seed ) != CURAND_STATUS_SUCCESS ) {
        assert(false);
    }

    uniform_fill( dGen, key_vec, N );
    uniform_fill( dGen, val_vec, N );

    periodic_sort( key_vec, val_vec, sorted_keys, sorted_vals );

    boost::property_tree::ptree log;
    bool valid = validate_sorted(key_vec, val_vec, sorted_keys, sorted_vals, log );

    if( !valid ) {
        std::cout << "FAILURE: Invalid Sorting" << std::endl;
        boost::property_tree::ptree input_log, output_log;

        dump( key_vec, val_vec, input_log );
        dump( sorted_keys, sorted_vals, output_log );

        log.add_child( "input", input_log);
        log.add_child( "output", output_log );

        log.put("seed", seed );
        log.put("N", N );
        log.put("Valid", valid );

        boost::property_tree::write_json( std::cerr, log );
    } else {
        std::cout << "SUCCESS" << std::endl;
    }

    return 0;
}
