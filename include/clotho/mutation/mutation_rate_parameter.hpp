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
#ifndef MUTATION_RATE_PARAMETER_HPP_
#define MUTATION_RATE_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

template < class RealType >
struct mutation_rate_parameter {
    RealType    m_mu;

    static constexpr RealType DEFAULT_MUTATION_RATE = 0.0001;

    mutation_rate_parameter( RealType m = DEFAULT_MUTATION_RATE ) :
        m_mu( m )
    { }

    mutation_rate_parameter( boost::property_tree::ptree & config ) :
        m_mu( DEFAULT_MUTATION_RATE )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( MUT_BLOCK_K, lconfig );

        m_mu = lconfig.get< RealType >( MU_K, m_mu );

        lconfig.put( MU_K, m_mu );
        config.put_child( MUT_BLOCK_K, lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( MU_K, m_mu );
        l.put_child( MUT_BLOCK_K, c );
    }

    virtual ~mutation_rate_parameter() {}
};

template < class RealType >
constexpr RealType mutation_rate_parameter< RealType >::DEFAULT_MUTATION_RATE;

#endif  // MUTATION_RATE_PARAMETER_HPP_
