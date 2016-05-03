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
#include "clotho/fitness/constant_fitness_metric.hpp"

const std::string CONSTANT_NAME = "constant";

constant_fitness_metric::constant_fitness_metric( real_type s ) {}

constant_fitness_metric::result_type constant_fitness_metric::operator()( real_type x ) {
    return 1.0;
}

constant_fitness_metric::result_type constant_fitness_metric::operator()( real_type x, real_type y ) {
    return 1.0;
}

const std::string constant_fitness_metric::name() const {
    return CONSTANT_NAME;
}

void constant_fitness_metric::log( std::ostream & out ) const {
    out << "{" << CONSTANT_NAME << "}\n";
}

constant_fitness_metric::~constant_fitness_metric() {}
