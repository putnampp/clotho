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
#ifndef CLOTHO_WEIGHT_DISTRIBUTION_HELPER_HPP_
#define CLOTHO_WEIGHT_DISTRIBUTION_HELPER_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/random/normal_distribution.hpp>
#include "clotho/random/normal_distribution_parameter.hpp"

namespace clotho {
namespace genetics {

template < class T >
struct weight_distribution_helper;


template < >
struct weight_distribution_helper< double > {
    typedef boost::random::normal_distribution< double > type;
    typedef typename type::param_type                   result_type;

    static result_type  makeParameter( boost::property_tree::ptree & config ) {
        normal_distribution_parameter< double > param( config );
        return result_type( param.m_mean, param.m_sigma );
    }
};

template < >
struct weight_distribution_helper< float > {
    typedef boost::random::normal_distribution< float >     type;
    typedef typename type::param_type                       result_type;

    static result_type  makeParameter( boost::property_tree::ptree & config ) {
        normal_distribution_parameter< float > param( config );
        return result_type( param.m_mean, param.m_sigma );
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_WEIGHT_DISTRIBUTION_HELPER_HPP_

