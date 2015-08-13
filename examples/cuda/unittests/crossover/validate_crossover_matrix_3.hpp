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
#ifndef VALIDATE_CROSSOVER_MATRIX_3_HPP_
#define VALIDATE_CROSSOVER_MATRIX_3_HPP_

#include "crossover_test.hpp"

#include "clotho/cuda/crossover/crossover_matrix_3.cuh"
#include "clotho/utility/log_helper.hpp"

bool validate( crossover_test< crossover< 3 > > &ct, boost::property_tree::ptree & err ) {

    typedef crossover< 3 > crossover_type;
    typedef crossover_test< crossover_type > test_type;

    typedef typename test_type::random_vector random_vector;
    
    typedef typename test_type::sequence_iterator sequence_iterator;
    sequence_iterator s_it = ct.sequences.begin();

    boost::property_tree::ptree seqs;
    int sidx = 0;
    while( s_it != ct.sequences.end() ) {
        boost::property_tree::ptree s;
        clotho::utility::add_value_array(s, s_it, s_it + ct.allele_list.size() / 32 );

        s_it += ct.allele_list.size() / 32;
    
        std::ostringstream oss;
        oss << sidx++;
        seqs.add_child( oss.str(), s);
    }

    err.add_child( "sequences", seqs );
    return false;
}

#endif  // VALIDATE_CROSSOVER_MATRIX_3_HPP_
