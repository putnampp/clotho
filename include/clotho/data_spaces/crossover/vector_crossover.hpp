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
#ifndef CLOTHO_VECTOR_CROSSOVER_HPP_
#define CLOTHO_VECTOR_CROSSOVER_HPP_

#include "clotho/data_spaces/crossover/block_crossover.hpp"
#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class ClassifierType, class PopulationType >
class vector_crossover : public block_crossover< Classifier, typename PopulationType::block_type >, public task {
public:
    typedef block_crossover< Classifier, typename PopulationType::block_type > base_type;

    typedef typename PopulationType::sequence_type  sequence_type;  // shared_ptr
    typedef typename PopulationType::block_type     block_type;

    typedef ClassifierType classifier_type;

    vector_crossover( const classifier_type & cfier, sequence_type top, sequence_type bottom, sequence_type res )

    void operator()() {
        typedef typename PopulationType::genetic_type::const_sequence_iterator const_iterator;

        const_iterator tb = top->begin_sequence(), te = top->end_sequence(), bb = bottom->begin_sequence(), be = bottom->end_sequence();

        unsigned int i = 0;
        while( true ) {
            if( tb == te ) {
                while( bb != be ) {
                    block_type o = this->crossover( typename base_type::bit_helper_type::ALL_UNSET, *bb++, i++ );
                    res->first.push_back(o);
                }
                break;

            } else if( bb == be ) {
                while( tb != te ) {
                    block_type o = this->crossover( *tb++, typename base_type::bit_helper_type::ALL_UNSET, i++ );
                    res->first.push_back(o);
                }
                break;
            }

            block_type t = *tb++;
            block_type b = *bb++;

            block_type o = this->crossover(t, b, i );

            res->first.push_back( o );
            ++i;
        }
    }

    virtual ~vector_crossover() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_VECTOR_CROSSOVER_HPP_

