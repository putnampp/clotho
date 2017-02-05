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
#ifndef CLOTHO_BLOCK_CROSSOVER_METHOD_HPP_
#define CLOTHO_BLOCK_CROSSOVER_METHOD_HPP_

#include "clotho/data_spaces/crossover/block_crossover.hpp"

namespace clotho {
namespace genetics {

template < class ClassifierType, class BlockType >
class block_crossover_method : public block_crossover< ClassifierType, BlockType > {
public:
    typedef block_crossover_method< ClassifierType, BlockType > self_type;

    typedef ClassifierType  classifier_type;
    typedef BlockType       block_type;
    typedef block_type *     block_pointer;

    typedef block_crossover< ClassifierType, BlockType >    operator_type;
    typedef typename operator_type::bit_helper_type         bit_helper_type;

    block_crossover_method( const classifier_type & cls ) : operator_type( cls ) {}

    void operator()( block_pointer p0_start, block_pointer p0_end, block_pointer p1_start, block_pointer p1_end, block_pointer o_start, block_pointer o_end ) {

        block_pointer offspring = o_start;
        unsigned int i = 0;
        while( true ) {
            if( p0_start == p0_end ) {
                while( p1_start != p1_end ) {
                    const block_type t = bit_helper_type::ALL_UNSET;
                    const block_type b = *p1_start++;
                    *offspring++ = this->crossover( t, b, i );

                    i += bit_helper_type::BITS_PER_BLOCK;
                }
                break;
            } else if( p1_start == p1_end ) {
                while( p0_start != p0_end ) {
                    const block_type t = *p0_start++;
                    const block_type b = bit_helper_type::ALL_UNSET;
                    *offspring++ = this->crossover( t, b, i );
                    i += bit_helper_type::BITS_PER_BLOCK;
                }
                break;
            }

            const block_type t = *p0_start++;
            const block_type b = *p1_start++;

            *offspring++ = this->crossover( t, b, i );
            i += bit_helper_type::BITS_PER_BLOCK;
        }

        while( offspring != o_end ) {
            *offspring++ = bit_helper_type::ALL_UNSET;
        }
    }

protected:
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_BLOCK_CROSSOVER_METHOD_HPP_
