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
#ifndef ELEMENT_CLASSIFIER_HPP_
#define ELEMENT_CLASSIFIER_HPP_

namespace clotho {
namespace classifier {


template < class Result = bool >
class element_classifier {
public:
    typedef Result  result_type;

    template < class Element >
    result_type operator()( const element_type & elem ) const {
        return false;
    }

    template < class ElementIterator >
    result_type operator()( ElementIterator elem_it, size_t idx) const {
        return operator()( *(elem_it + idx) );
    }
};

}   // namespace classifier {
}   // namespace clotho {

#endif  // ELEMENT_CLASSIFIER_HPP_
