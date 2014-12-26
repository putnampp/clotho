#ifndef NORMAL_FITNESS_METRIC_HPP_
#define NORMAL_FITNESS_METRIC_HPP_

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <cmath>

/**
 * univariate normal distribution function
 *
 * Note: for k-dimensional vector the fitness is that of the first component
 */
class normal_fitness_metric {
public:
    typedef double      real_type;
    typedef real_type   result_type;

    normal_fitness_metric( real_type mu = 0., real_type sigma = 1. ) :
        m_mean( mu )
        , m_sigma( sigma )
        , _coeff( (1.0/sigma) * boost::math::double_constants::one_div_root_two_pi)
        , _denom( 2.0 * sigma * sigma )
    {}

    result_type operator()( real_type x ) {
        double ex = (x - m_mean);
        ex *= ex;
        ex /= _denom;
        return _coeff * exp( -ex );
    }

    result_type operator()( real_type x, real_type mu, real_type sigma ) {
        double ex = (x - mu);
        ex *= ex;
        ex /= (2.0 * sigma * sigma );

        double coeff = (1.0/sigma) * boost::math::double_constants::one_div_root_two_pi;

        return coeff * exp(-ex);
    }

    result_type operator()( const std::vector< real_type > & multi_variate ) {
        if( multi_variate.empty() ) return operator()( 0. );

        return operator()( multi_variate.front() );
    }

    result_type operator()( const std::vector< real_type > & multi_variate, real_type mu, real_type sigma ) {
        if( multi_variate.empty() ) return operator()( 0., mu, sigma );

        return operator()( multi_variate.front(), mu, sigma );
    }

protected:
    real_type m_mean;
    real_type m_sigma;

    real_type _coeff, _denom;
};

#endif  // NORMAL_FITNESS_METRIC_HPP_
