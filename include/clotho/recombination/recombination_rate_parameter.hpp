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
#ifndef RECOMBINATION_RATE_PARAMETER_HPP_
#define RECOMBINATION_RATE_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

template < class RealType >
struct recombination_rate_parameter {
    RealType    m_rho;

    static constexpr RealType DEFAULT_RECOMB_RATE = 0.0001;

    recombination_rate_parameter( RealType r = DEFAULT_RECOMB_RATE ) :
        m_rho( r )
    {}

    recombination_rate_parameter( boost::property_tree::ptree & config ) :
        m_rho( DEFAULT_RECOMB_RATE )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( REC_BLOCK_K, lconfig );

        m_rho = lconfig.get< RealType >( RHO_K, m_rho );

        lconfig.put( RHO_K, m_rho );
        config.put_child( REC_BLOCK_K, lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( RHO_K, m_rho );
        l.put_child( REC_BLOCK_K, c );
    }

    virtual ~recombination_rate_parameter() {}
};

template < class RealType >
constexpr RealType  recombination_rate_parameter< RealType >::DEFAULT_RECOMB_RATE;

#endif  // RECOMBINATION_RATE_PARAMETER_HPP_
