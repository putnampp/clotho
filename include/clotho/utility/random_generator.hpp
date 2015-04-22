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
#ifndef RANDOM_GENERATOR_HPP_
#define RANDOM_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

template < class URNG, class Element >
class random_generator {
public:
    typedef random_generator< URNG, Element > self_type;
    typedef Element                           result_type;

    //random_generator( URNG & rng, boost::property_tree::ptree & config ) : m_rng( &rng ) {}

    result_type operator()() {
        return result_type( );
    }

protected:
    URNG * m_rng;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // RANDOM_GENERATOR_HPP_
