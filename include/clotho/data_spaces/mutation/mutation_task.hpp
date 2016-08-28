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
#ifndef CLOTHO_MUTATION_TASK_HPP_
#define CLOTHO_MUTATION_TASK_HPP_

namespace clotho {
namespace genetics {

template < class RNG >
class mutation_task : public task {
public:

    typedef RNG random_generator_type;

    template < class SeedSeq >
    mutation_task( SeedSeq & seed, double mu ) :
        m_rand( seed )
        , m_allele_gen( &m_rand, mu )
    {}

    void operator()() {
    }

    virtual ~mutation_task() {}

protected:

    random_generator_type   m_rand;
};

}   // namespace genetics
}   // namespace clotho
#endif  //CLOTHO_MUTATION_TASK_HPP_
