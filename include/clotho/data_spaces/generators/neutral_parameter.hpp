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
#ifndef CLOTHO_NEUTRAL_PARAMETER_HPP_
#define CLOTHO_NEUTRAL_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

namespace clotho {
namespace genetics {

template < class RealType >
struct neutral_parameter {
    RealType    m_p;

    static constexpr RealType DEFAULT_NEUTRALITY = 0.5;

    neutral_parameter( RealType m = DEFAULT_NEUTRALITY ) :
        m_p( m )
    { }

    neutral_parameter( boost::property_tree::ptree & config ) :
        m_p( DEFAULT_NEUTRALITY )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( NEUTRAL_BLOCK_K, lconfig );

        m_p = lconfig.get< RealType >( P_K, m_p );

        lconfig.put( P_K, m_p );
        config.put_child( NEUTRAL_BLOCK_K, lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( P_K, m_p );
        l.put_child( NEUTRAL_BLOCK_K, c );
    }

    virtual ~neutral_parameter() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_NEUTRAL_PARAMETER_HPP_
