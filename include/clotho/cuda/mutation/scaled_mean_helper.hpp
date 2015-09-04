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
#ifndef SCALED_MEAN_HELPER_HPP_
#define SCALED_MEAN_HELPER_HPP_

template < class RealType, class OrderTag >
struct scaled_mean_helper {
    static RealType get( RealType m ) {
        return m;
    }
};

#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

template < class RealType, class IntType >
struct scaled_mean_helper< RealType, unit_ordered_tag< IntType > > {

    static RealType get( RealType m ) {
        RealType s = (m / ((RealType) unit_ordered_tag< IntType >::OBJECTS_PER_INT ));
        return s;
    }
};

#endif  // SCALED_MEAN_HELPER_HPP_
