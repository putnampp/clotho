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
#ifndef CLOTHO_PAIRWISE_DIFFERENCE_DEF_HPP_
#define CLOTHO_PAIRWISE_DIFFERENCE_DEF_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace genetics {

template < class SequenceSpaceType >
class pairwise_difference {
public:
    typedef SequenceSpaceType   space_type;

    void evaluate( space_type & ss, boost::property_tree::ptree & results );

    virtual ~pairwise_difference() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PAIRWISE_DIFFERENCE_DEF_HPP_
