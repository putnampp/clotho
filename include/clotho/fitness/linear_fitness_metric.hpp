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
#ifndef LINEAR_FITNESS_METRIC_HPP_
#define LINEAR_FITNESS_METRIC_HPP_

#include <vector>
#include "clotho/genetics/ifitness.hpp"

extern const std::string LINEAR_FIT_NAME;

class linear_fitness_metric : public ifitness {
public:
    typedef double      real_type;
    typedef real_type   result_type;

    linear_fitness_metric( real_type a = 1., real_type b = 0. ) :
        m_val( c ) {
    }

    result_type operator()() {
        return m_B;
    }

    result_type operator()( real_type x ) {
        return m_A * x + m_B;
    }

    result_type operator()( const std::vector< real_type > & multi_variate ) {
        return m_B;
    }

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~linear_fitness_metric();
protected:
    real_type m_A, m_B;
};

#endif  // LINEAR_FITNESS_METRIC_HPP_
