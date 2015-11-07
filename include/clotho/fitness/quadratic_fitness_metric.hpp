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
#ifndef QUADRATIC_FITNESS_METRIC_HPP_
#define QUADRATIC_FITNESS_METRIC_HPP_

#include "clotho/fitness/ifitness.hpp"
#include "clotho/utility/clotho_strings.hpp"

/**
 * Fitness is the scaled phenotype
 *
 * f( x ) = g( x; s ), if g(x ; s) > 0;
 *        = 0, otherwise
 *
 * g( x; s ) = 1.0 - (x / s)^2
 *
 * x - specific trait's phenotype   (double)
 * s - scaling factor               (double)
 *
 */
class quadratic_fitness_metric : public ifitness {
public:
    typedef double real_type;
    typedef real_type result_type;

    quadratic_fitness_metric( real_type s = 1. );

    result_type operator()( real_type x );
    result_type operator()( real_type x, real_type s );

    inline result_type operator()( const std::vector< real_type > & multi_variate ) {
        return ((multi_variate.empty()) ? 0.0 :  operator()( multi_variate.front() ));
    }

    inline result_type operator()( const std::vector< real_type > & multi_variate, real_type s ) {
        return ((multi_variate.empty()) ? 0.0 : operator()( multi_variate.front(), s ));
    }

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~quadratic_fitness_metric();

protected:
    real_type m_std;
};

#endif  // QUADRATIC_FITNESS_METRIC_HPP_
