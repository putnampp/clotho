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
#ifndef NORMAL_FITNESS_METRIC_HPP_
#define NORMAL_FITNESS_METRIC_HPP_

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <cmath>

#include "clotho/fitness/ifitness.hpp"

extern const std::string NORM_NAME;

/**
 * univariate normal distribution function
 *
 * Note: for k-dimensional trait vector the fitness is based on the phenotype of
 * only the first trait
 */
class normal_fitness_metric : public ifitness {
public:
    typedef double      real_type;
    typedef real_type   result_type;

    normal_fitness_metric( real_type mu = 0., real_type sigma = 1. );

    result_type operator()( real_type x );
    result_type operator()( real_type x, real_type mu, real_type sigma );

    inline result_type operator()( const std::vector< real_type > & multi_variate ) {
        return ((multi_variate.empty()) ? operator()( 0. ) : operator()( multi_variate.front() ));
    }

    inline result_type operator()( const std::vector< real_type > & multi_variate, real_type mu, real_type sigma ) {
        return ((multi_variate.empty()) ? operator()( 0., mu, sigma ) : operator()( multi_variate.front(), mu, sigma ));
    }

    inline result_type operator()( real_type * first, real_type * last ) {
        return (first == last) ? operator()( 0. ) : operator()( *first );
    }

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~normal_fitness_metric();

protected:
    real_type m_mean;
    real_type m_sigma;

    real_type _coeff, _denom;
};

#endif  // NORMAL_FITNESS_METRIC_HPP_
