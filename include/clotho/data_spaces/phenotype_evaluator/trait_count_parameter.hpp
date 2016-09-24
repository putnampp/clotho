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
#ifndef TRAIT_COUNT_PARAMETER_HPP_
#define TRAIT_COUNT_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

struct trait_count_parameter {
    unsigned int m_trait_count;

    static const unsigned int DEFAULT_TRAIT_COUNT = 1;

    trait_count_parameter( unsigned int N = DEFAULT_TRAIT_COUNT ) :
        m_trait_count(N)
    {}

    trait_count_parameter( boost::property_tree::ptree & config ) :
        m_trait_count( DEFAULT_TRAIT_COUNT )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( "traits", lconfig );

        m_trait_count = lconfig.get< unsigned int >( "N", m_trait_count );

        lconfig.put( "N", m_trait_count );
        config.put_child( "traits", lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( "N", m_trait_count );
        l.put_child( "traits", c );
    }

    virtual ~trait_count_parameter() {}
};

#endif  // TRAIT_COUNT_PARAMETER_HPP_
