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
#ifndef CROSSOVER_MATRIX_DEF_HPP_
#define CROSSOVER_MATRIX_DEF_HPP_

#include "clotho/cuda/compute_capability.hpp"

template < unsigned int Version >
class crossover {
public:
    typedef double                      real_type;
    typedef double                      allele_type;
    typedef unsigned int                event_count_type;
    typedef unsigned int                int_type;
    typedef unsigned int                size_type;
    typedef compute_capability< 3, 0 >  comp_cap_type;

    static const unsigned int           MAX_EVENTS;
    
    void operator()(  real_type         * rand_pool
                    , allele_type       * allele_list
                    , event_count_type  * event_list
                    , int_type          * sequences
                    , size_type nSequences
                    , size_type nAlleles
                    , size_type sequence_width );

    void get_state( boost::property_tree::ptree & s ) {
        s.put("crossover.version", Version );
    }

    virtual ~crossover() {}
};

#endif  // CROSSOVER_MATRIX_DEF_HPP_
