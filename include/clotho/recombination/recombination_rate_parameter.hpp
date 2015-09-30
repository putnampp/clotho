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

template < class RealType >
struct recombination_rate_parameter {
    RealType    m_rho;

    recombination_rate_parameter( boost::property_tree::ptree & config ) :
        m_rho( 0.001 )
    {
        boost::property_tree::ptree lconfig;

        if( config.get_child_optional( "recombination" ) != boost::none ) {
            lconfig = config.get_child( "recombination" );
        }

        if( lconfig.get_child_optional( "rho" ) == boost::none ) {
            lconfig.put("rho", m_rho );
        } else {
            m_rho = lconfig.get< RealType >( "rho", m_rho );
        }

        config.put_child( "recombination", lconfig );
    }

    virtual ~recombination_rate_parameter() {}
};

#endif  // RECOMBINATION_RATE_PARAMETER_HPP_
