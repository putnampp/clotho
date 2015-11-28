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
#ifndef SEED_PARAMETER_HPP_
#define SEED_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"
#include "clotho/utility/timer.hpp"

template < class SeedType = unsigned long long >
struct seed_parameter {
    typedef SeedType seed_type;

    static const seed_type DEFAULT_SEED = 0;

    seed_type m_seed;

    seed_parameter( seed_type s ) :
        m_seed( s )
    {}

    seed_parameter( boost::property_tree::ptree & config ) :
        m_seed( DEFAULT_SEED )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( RNG_BLOCK_K, lconfig );

        m_seed = lconfig.get< seed_type >( SEED_K, m_seed );

        if( m_seed == DEFAULT_SEED ) {
            m_seed = clotho::utility::clock_type::now().time_since_epoch().count();
        }

        lconfig.put( SEED_K, m_seed );
        config.put_child( RNG_BLOCK_K, lconfig );
    }

    seed_parameter( const seed_parameter< seed_type > & rhs ) : m_seed( rhs.m_seed ) {}

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( SEED_K, m_seed );
        l.put_child( RNG_BLOCK_K, c );
    }

    seed_parameter< seed_type > & operator=( const seed_parameter< seed_type > & rhs ) {
        m_seed = rhs.m_seed;
        return *this;
    }

    virtual ~seed_parameter() {}
};

template < class SeedType >
const SeedType seed_parameter< SeedType >::DEFAULT_SEED;

#endif  // SEED_PARAMETER_HPP_
