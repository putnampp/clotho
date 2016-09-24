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
#ifndef CLOTHO_PAIRWISE_PHENOTYPE_REDUCER_HPP_
#define CLOTHO_PAIRWISE_PHENOTYPE_REDUCER_HPP_

namespace clotho {
namespace genetics {

struct pairwise_phenotype_reducer {
    template < class PhenotypeType, class ResultType >
    void operator()( ResultType * res, PhenotypeType * s0, PhenotypeType * s1, unsigned int N ) {
        unsigned int i = 0;
        while( i < N ) {
            res[i] = s0[i] + s1[i];
            ++i;
        }
    }
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_PAIRWISE_PHENOTYPE_REDUCER_HPP_
