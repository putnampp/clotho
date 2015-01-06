#ifndef NORMAL_FITNESS_METRIC_HPP_
#define NORMAL_FITNESS_METRIC_HPP_

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <cmath>

#include "clotho/genetics/ifitness.hpp"

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

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~normal_fitness_metric();

protected:
    real_type m_mean;
    real_type m_sigma;

    real_type _coeff, _denom;
};

#endif  // NORMAL_FITNESS_METRIC_HPP_
