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
#ifndef FITNESS_METHOD_HPP_
#define FITNESS_METHOD_HPP_

namespace clotho {
namespace fitness {

struct multiplicative_heterozygous_tag {};
struct multiplicative_homozygous_tag {};

struct additive_heterozygous_tag {};
struct additive_homozygous_tag {};

template < class Result, class Tag >
class fitness_method {
public:
    typedef Result result_type;

    template < class Element >
    void operator()( result_type & res, const Element & elem, result_type scale = 1. ) {}
};

}   // namespace fitness {
}   // namespace clotho {

#endif  // FITNESS_METHOD_HPP_
