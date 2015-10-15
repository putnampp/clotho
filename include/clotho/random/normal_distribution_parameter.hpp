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
#ifndef NORMAL_DISTRIBUTION_PARAMETER_HPP_
#define NORMAL_DISTRIBUTION_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"


template < class RealType >
struct normal_distribution_parameter {
    static constexpr RealType   DEFAULT_MEAN = 0.0;
    static constexpr RealType   DEFAULT_SIGMA = 1.0;

    RealType m_mean, m_sigma;

    normal_distribution_parameter( RealType m = DEFAULT_MEAN, RealType s = DEFAULT_SIGMA ):
        m_mean(m)
        , m_sigma(s)
    {}

    normal_distribution_parameter( boost::property_tree::ptree & config ) :
        m_mean( DEFAULT_MEAN )
        , m_sigma( DEFAULT_SIGMA )
    {
        m_mean = config.get< RealType >( MEAN_K, m_mean );
        m_sigma = config.get< RealType >( SIGMA_K, m_sigma );

        config.put( MEAN_K, m_mean );
        config.put( SIGMA_K, m_sigma );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        l.put( MEAN_K, m_mean );
        l.put( SIGMA_K, m_sigma );
    }

    virtual ~normal_distribution_parameter() {}
};

#endif  // NORMAL_DISTRIBUTION_PARAMETER_HPP_
