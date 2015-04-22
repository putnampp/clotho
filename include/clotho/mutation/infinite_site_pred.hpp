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
#ifndef CLOTHO_INFINITE_SITE_PRED_HPP_
#define CLOTHO_INFINITE_SITE_PRED_HPP_

#include "clotho/mutation/no_duplicate_pred.hpp"

namespace clotho {
namespace mutations {

template < class Set >
struct infinite_site_pred {
    typedef no_duplicate_pred< Set > type;
};

}   // namespace mutations
}   // namespace clotho

#endif  // CLOTHO_INFINITE_SITE_PRED_HPP_
