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
#ifndef IS_EVEN_TAG_HPP_
#define IS_EVEN_TAG_HPP_

namespace clotho {
namespace classifiers {
namespace tags {

struct is_even_tag {
    template < class IntType>
    static bool eval( IntType t ) {
        return !(t & 1);
    }
};

}   // namespace tags {
}   // namespace classifiers {
}   // namespace clotho {

#endif  // IS_EVEN_TAG_HPP_
