#ifndef CONSTANT_FITNESS_METRIC_HPP_
#define CONSTANT_FITNESS_METRIC_HPP_

class constant_fitness_metric {
public:
    typedef double      real_type;
    typedef real_type   result_type;

    constant_fitness_metric( real_type c = 1. ) :
        m_val( c ) {
    }

    result_type operator()() {
        return m_val;
    }

    result_type operator()( real_type x ) {
        return m_val;
    }

    result_type operator()( const std::vector< real_type > & multi_variate ) {
        return m_val;
    }

protected:
    real_type m_val;
};

#endif  // CONSTANT_FITNESS_METRIC_HPP_
