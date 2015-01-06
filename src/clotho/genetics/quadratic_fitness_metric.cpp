#include "clotho/genetics/quadratic_fitness_metric.hpp"
#include <cassert>

const std::string QUAD_NAME = "quadratic";

quadratic_fitness_metric::quadratic_fitness_metric( real_type s ) :
    m_scale(s)
{
    assert( m_scale != 0 );
}

quadratic_fitness_metric::result_type quadratic_fitness_metric::operator()( real_type x ) {
    result_type res = x / m_scale;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res));
}

quadratic_fitness_metric::result_type quadratic_fitness_metric::operator()( real_type x, real_type s ) {
    assert( s != 0 );

    result_type res = x / s;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res ));
}

const std::string quadratic_fitness_metric::name() const {
    return QUAD_NAME;
}

void quadratic_fitness_metric::log( std::ostream & out ) const {
    out << "{" << QUAD_NAME << "," << m_scale << "}\n";
}

quadratic_fitness_metric::~quadratic_fitness_metric() {}
