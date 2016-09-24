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
#ifndef SEQUENCE_BIAS_PARAMETER_HPP_
#define SEQUENCE_BIAS_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

template < class RealType >
struct sequence_bias_parameter {
    RealType m_bias;
    static constexpr RealType DEFAULT_BIAS = 0.5;

    sequence_bias_parameter( RealType b = DEFAULT_BIAS ) : 
        m_bias( b ) 
    {}

    sequence_bias_parameter( boost::property_tree::ptree & config ) : 
        m_bias( DEFAULT_BIAS ) 
    {

        boost::property_tree::ptree lconfig, bblock;

        lconfig = config.get_child( REC_BLOCK_K, lconfig );
        bblock = lconfig.get_child( SEQUENCE_BIAS_BLOCK_K, bblock );

        m_bias = bblock.get< RealType >( P_K, m_bias);

        bblock.put( P_K, m_bias );
        lconfig.put_child( SEQUENCE_BIAS_BLOCK_K, bblock);
        config.put_child( REC_BLOCK_K, lconfig );
    }

    virtual ~sequence_bias_parameter() {}
};

#endif  // SEQUENCE_BIAS_PARAMETER_HPP_
