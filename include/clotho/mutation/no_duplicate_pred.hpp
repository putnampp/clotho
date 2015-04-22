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
#ifndef CLOTHO_NO_DUPLICATE_HPP_
#define CLOTHO_NO_DUPLICATE_HPP_

namespace clotho {
namespace mutations {

template < class Set >
class no_duplicate_pred {
public:
    typedef Set                 set_type;
    typedef typename Set::element_type   element_type;

    no_duplicate_pred( set_type * elements ) :
        m_elements( elements ) {
    }

    bool operator()( const element_type & elem ) {
        assert(false );
    }

    virtual ~no_duplicate_pred() {}
protected:
    set_type    * m_elements;
};

}   // namespace mutations {
}   // namespace clotho {

#endif  // CLOTHO_NO_DUPLICATE_HPP_
