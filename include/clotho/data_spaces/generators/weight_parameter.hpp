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
#ifndef CLOTHO_WEIGHT_PARAMETER_HPP_
#define CLOTHO_WEIGHT_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/utility/clotho_strings.hpp"

#include "clotho/data_spaces/generators/weight_distribution_helper.hpp"


namespace clotho {
namespace genetics {

template < class WeightType >
struct weight_parameter {
    typedef WeightType  weight_type;
    typedef  weight_distribution_helper< weight_type > dist_helper_type;
    typedef typename dist_helper_type::type distribution_type;

    distribution_type   m_dist;

    weight_parameter( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( TRAIT_BLOCK_K, lconfig );
        m_dist.param( dist_helper_type::makeParameter( lconfig) );

        config.put_child( TRAIT_BLOCK_K, lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        dist_helper_type::makeParameter( c );
        l.put_child( TRAIT_BLOCK_K, c );
    }

    virtual ~weight_parameter() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_WEIGHT_PARAMETER_HPP_
