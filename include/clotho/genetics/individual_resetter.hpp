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
#ifndef INDIVIDUAL_RESETTER_HPP_
#define INDIVIDUAL_RESETTER_HPP_

template < class IndividualType >
struct individual_resetter;

#include "resetter.hpp"

template < class SequenceType >
struct individual_resetter< std::pair< SequenceType, SequenceType > > {
    typedef std::pair< SequenceType, SequenceType > individual_type;
    typedef resetter< SequenceType > resetter_type;

    void operator()( individual_type & ind ) {
        resetter_type::reset( ind.first );
        resetter_type::reset( ind.second  );;
    }
};

#endif  // INDIVIDUAL_RESETTER_HPP_
