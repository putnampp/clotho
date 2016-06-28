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
#ifndef CLOTHO_LINEAR_COMBINATION_METHOD_HPP_
#define CLOTHO_LINEAR_COMBINATION_METHOD_HPP_

#include <vector>
#include "clotho/data_spaces/accumulator_helper_of.hpp"

namespace clotho {
namespace genetics {

template < class TraitAccumType, class ResultType >
struct linear_combination {

    typedef TraitAccumType                                  trait_accum_type;
    typedef typename trait_accum_type::trait_vector_type    trait_type;
    typedef typename trait_accum_type::weight_pointer       weight_pointer;
    typedef ResultType                                      result_type;

    result_type operator()( std::shared_ptr< trait_type > s0, std::shared_ptr< trait_type > s1 ) {
        typedef typename trait_type::iterator iterator;

        iterator b = s0->begin(), e = s0->end();
        result_type res = 0.0;

        while( b != e ) {
            res += *b++;
        }

        b = s1->begin(); e = s1->end();

        while( b != e ) {
            res += *b++;
        }

        return res;
    }

    void operator()( weight_pointer a_first, weight_pointer a_last, weight_pointer b_first, weight_pointer b_last, weight_pointer res ) {
#ifdef DEBUGGING
        std::cerr << "Linear combination method" << std::endl;
        std::cerr << "Trait count A: " << (a_last - a_first) << std::endl;
        std::cerr << "Trait count B: " << (b_last - b_first) << std::endl;
#endif  // DEBUGGING
        int i = 0, M = (a_last - a_first);

        assert( M == (b_last - b_first) );
        while( i < M ) {
            res[ i ] = a_first[ i ] + b_first[ i ];
            ++i;
        }
    }
};

template < class TraitAccumType, class RealType >
struct linear_combination< TraitAccumType, std::vector< RealType > > {

    typedef TraitAccumType                                  trait_accum_type;
    typedef typename trait_accum_type::trait_vector_type    trait_type;
    typedef std::vector< RealType >                         result_type;

    result_type operator()( std::shared_ptr < trait_type >  s0, std::shared_ptr <trait_type>  s1 ) {
        typedef typename trait_type::iterator iterator;

        iterator b = s0->begin(), e = s0->end();
        result_type res;

        while( b != e ) {
            res.push_back( *b++ );
        }

        b = s1->begin(); e = s1->end();

        size_t i = 0;
        while( b != e ) {
            res[ i++ ] += *b++;
        }

        return res;
    }
};
template < class TraitAccumType, class ResultType >
struct accumulator_helper_of< linear_combination< TraitAccumType, ResultType > > {
    typedef TraitAccumType          type;
    typedef ResultType              result_type;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_LINEAR_COMBINATION_METHOD_HPP_
