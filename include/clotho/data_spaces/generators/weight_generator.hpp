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
#ifndef CLOTHO_WEIGHT_GENERATOR_HPP_
#define CLOTHO_WEIGHT_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/generators/weight_parameter.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class WeightType >
class weight_generator : public weight_parameter< WeightType > {
public:
    typedef weight_parameter< WeightType >      base_type;

    typedef RNG                                 random_engine_type;
    typedef typename base_type::weight_type      weight_type;


    weight_generator( random_engine_type * rng, boost::property_tree::ptree & config ) : 
        base_type ( config )
        , m_rand( rng )
    { }

    weight_type   operator()() {
        return this->m_dist( * m_rand );
    }

    virtual ~weight_generator() {}

protected:
    random_engine_type  * m_rand;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_WEIGHT_GENERATOR_HPP_
