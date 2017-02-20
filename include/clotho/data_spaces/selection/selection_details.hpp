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
#ifndef CLOTHO_SELECTION_DETAILS_HPP_
#define CLOTHO_SELECTION_DETAILS_HPP_

namespace clotho {
namespace genetics {

template < class RNG, class IDType = unsigned int >
struct selection_details {
    typedef RNG                                     random_engine_type;
    typedef IDType                                  individual_id_type;
    typedef std::pair< individual_id_type, individual_id_type > parent_pair;

    random_engine_type * m_rand;

    selection_details( random_engine_type * rng ) : m_rand(rng) {}

    virtual parent_pair operator()() = 0;

    virtual ~selection_details() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_SELECTION_DETAILS_HPP_
