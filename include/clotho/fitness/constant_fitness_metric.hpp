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
#ifndef CONSTANT_FITNESS_METRIC_HPP_
#define CONSTANT_FITNESS_METRIC_HPP_

#include "clotho/fitness/ifitness.hpp"

extern const std::string CONSTANT_NAME;

class constant_fitness_metric : public ifitness {
public:

    typedef double real_type;
    typedef real_type result_type;

    constant_fitness_metric( real_type s );

    result_type operator()( real_type x );
    result_type operator()( real_type x, real_type s );

    inline result_type operator()( const std::vector< real_type > & multi_variate ) {
        return 1.0;
    }

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~constant_fitness_metric();
};

#endif  // CONSTANT_FITNESS_METRIC_HPP_
