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
#ifndef SEQUENCE_HELPER_HPP_
#define SEQUENCE_HELPER_HPP_

#include "clotho/utility/iterator_helper.hpp"

template < class IndividualType >
struct sequence_helper {
    typedef typename IndividualType::sequence_type type;
    typedef clotho::utility::iterator_helper< IndividualType > iterator_helper_type;
};

template < class SequenceType >
struct sequence_helper< std::pair< SequenceType, SequenceType > > {
    typedef SequenceType        sequence_type;
    typedef clotho::utility::iterator_helper< std::pair< SequenceType, SequenceType > > iterator_helper_type;
};

#endif  // SEQUENCE_HELPER_HPP_
