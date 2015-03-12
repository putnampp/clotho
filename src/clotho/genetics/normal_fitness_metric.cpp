#include "clotho/genetics/normal_fitness_metric.hpp"

const std::string NORM_NAME = "gaussian";

normal_fitness_metric::normal_fitness_metric( real_type mu, real_type sigma ) :
    m_mean( mu )
    , m_sigma( sigma )
    , _coeff( (1.0/sigma) * boost::math::double_constants::one_div_root_two_pi)
    , _denom( 2.0 * sigma * sigma ) {
}

normal_fitness_metric::result_type normal_fitness_metric::operator()( real_type x ) {
    double ex = (x - m_mean);
    ex *= ex;
    ex /= (_denom);

    return _coeff * exp( -ex );
}

normal_fitness_metric::result_type normal_fitness_metric::operator()( real_type x, real_type mu, real_type sigma ) {
    double ex = (x - mu);
    ex *= ex;
    ex /= (2.0 * sigma * sigma );

    double coeff = (1.0/sigma) * boost::math::double_constants::one_div_root_two_pi;

    return coeff * exp(-ex);
}

const std::string normal_fitness_metric::name() const {
    return NORM_NAME;
}

void normal_fitness_metric::log( std::ostream & out ) const {
    out << "{" << NORM_NAME
        << ", " << m_mean
        << ", " << m_sigma
        << ", " << _coeff
        << ", " << _denom
        << "}\n";
}

normal_fitness_metric::~normal_fitness_metric() {}
