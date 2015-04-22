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
#ifndef INSPECT_METHODS_HPP_
#define INSPECT_METHODS_HPP_

namespace clotho {
namespace recombine {
namespace inspection {
namespace tag {

struct copy_matching_classify_mismatch {
    template < class Block >
    static Block select( const Block b0, const Block b1 ) {
        return (b0 ^ b1);
    }
};

struct classify_all {
    template < class Block >
    static Block select( const Block b0, const Block b1 ) {
        return (b0 | b1);
    }
};

}   // namespace tag
}   // namespace inspection
}   // namespace recombine
}   // namespace clotho

#endif  // INSPECT_METHODS_HPP_
