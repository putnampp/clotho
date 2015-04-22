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
#ifndef INDIVIDUAL_REPRODUCTION_DEF_HPP_
#define INDIVIDUAL_REPRODUCTION_DEF_HPP_

template < class IndividualType, class MutationGenerator, class RecombinationGenerator, class Tag >
class individual_reproduction {
public:
    typedef IndividualType individual_type;
    typedef IndividualType result_type;

    result_type operator()( individual_type & p0, individual_type & p1 ) {
        assert(false);
    }
};

#endif  // INDIVIDUAL_REPRODUCTION_DEF_HPP_
