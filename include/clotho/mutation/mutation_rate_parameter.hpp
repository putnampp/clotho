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

template < class RealType >
struct mutation_rate_parameter {
    RealType    m_mu;

    mutation_rate_parameter( boost::property_tree::ptree & config ) :
        m_mu( 0.001 )
    {
        boost::property_tree::ptree lconfig;

        if( config.get_child_optional( "mutation" ) != boost::none ) {
            lconfig = config.get_child( "mutation" );
        }

        if( lconfig.get_child_optional( "mutation_per_sequence" ) == boost::none ) {
            lconfig.put("mutation_per_sequence", m_mu );
        } else {
            m_mu = lconfig.get< RealType >( "mutation_per_sequence", m_mu );
        }

        config.put_child( "mutation", lconfig );
    }

    virtual ~mutation_rate_parameter() {}
};

#endif  // MUTATION_RATE_PARAMETER_HPP_
