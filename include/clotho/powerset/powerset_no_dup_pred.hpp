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
#ifndef CLOTHO_POWERSET_NO_DUP_PRED_HPP_
#define CLOTHO_POWERSET_NO_DUP_PRED_HPP_

#include "clotho/powerset/powerset.hpp"
#include "clotho/mutation/no_duplicate_pred.hpp"

namespace clotho {
namespace mutations {

template < class Element, class Subset, class Block, class BlockMap, class ElementKeyer >
class no_duplicate_pred< clotho::powersets::powerset< Element, Subset, Block, BlockMap, ElementKeyer> > {
public:
    typedef clotho::powersets::powerset< Element, Subset, Block, BlockMap, ElementKeyer> set_type;
    typedef Element   element_type;

    no_duplicate_pred( set_type * elements ) :
        m_elements( elements ) {
    }

    bool operator()( const element_type & elem ) {
        std::pair< typename set_type::element_index_type, bool > res = m_elements->find_or_create( elem );
        return res.second;
    }

    virtual ~no_duplicate_pred() {}
protected:
    set_type    * m_elements;
};

}   // namespace mutations {
}   // namespace clotho {

#endif  // CLOTHO_POWERSET_NO_DUP_PRED_HPP_
