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
#ifndef RANDOM_FITNESS_METRIC_HPP_
#define RANDOM_FITNESS_METRIC_HPP_

#include "clotho/fitness/ifitness.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

extern const std::string RAND_NAME;

class random_fitness_metric : public ifitness {
public:
    typedef double real_type;
    typedef real_type result_type;
    typedef unsigned int seed_type;

    typedef boost::random::mt19937      rng_type;
    typedef boost::random::uniform_01   distribution_type;

    random_fitness_metric( seed_type s = 0  );

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

    virtual ~random_fitness_metric();

protected:
    real_type m_scale;
};

#endif  // RANDOM_FITNESS_METRIC_HPP_
