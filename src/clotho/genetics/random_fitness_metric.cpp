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
#include "clotho/genetics/random_fitness_metric.hpp"
#include <cassert>

const std::string RAND_NAME = "random";

random_fitness_metric::random_fitness_metric( real_type s ) :
    m_scale(s) {
    assert( m_scale != 0 );
}

random_fitness_metric::result_type random_fitness_metric::operator()( real_type x ) {
    result_type res = x / m_scale;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res));
}

random_fitness_metric::result_type random_fitness_metric::operator()( real_type x, real_type s ) {
    assert( s != 0 );

    result_type res = x / s;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res ));
}

const std::string random_fitness_metric::name() const {
    return RAND_NAME;
}

void random_fitness_metric::log( std::ostream & out ) const {
    out << "{" << RAND_NAME << "," << m_scale << "}\n";
}

random_fitness_metric::~random_fitness_metric() {}
