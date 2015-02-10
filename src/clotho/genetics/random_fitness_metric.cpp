#include "clotho/genetics/random_fitness_metric.hpp"
#include <cassert>

const std::string RAND_NAME = "random";

random_fitness_metric::random_fitness_metric( real_type s ) :
    m_scale(s)
{
    assert( m_scale != 0 );
}

random_fitness_metric::result_type random_fitness_metric::operator()( real_type x ) {
    result_type res = x / m_scale;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res));
}

random_fitness_metric::result_type random_fitness_metric::operator()( real_type x, real_type s ) {
    assert( s != 0 );

    result_type res = x / s;
    res *= res;

    return (( res >= 1.0 ) ? 0.0 : (1.0 - res ));
}

const std::string random_fitness_metric::name() const {
    return RAND_NAME;
}

void random_fitness_metric::log( std::ostream & out ) const {
    out << "{" << RAND_NAME << "," << m_scale << "}\n";
}

random_fitness_metric::~random_fitness_metric() {}
