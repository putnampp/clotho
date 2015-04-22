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
#ifndef QTL_WEIGHT_ITERATOR_HPP_
#define QTL_WEIGHT_ITERATOR_HPP_

#include "clotho/utility/iterator_helper.hpp"
#include "qtl_allele.h"

namespace clotho {
namespace utility {

template < >
struct iterator_helper< qtl_allele > {
    typedef qtl_allele::weight_iterator type;

    static type make_first( qtl_allele & all ) {
        return all.begin();
    }

    static type make_last( qtl_allele & all ) {
        return all.end();
    }
};

}   // namespace utility {
}   // namespace clotho {
#endif  // QTL_WEIGHT_ITERATOR_HPP_
