#ifndef QUADRATIC_FITNESS_METRIC_HPP_
#define QUADRATIC_FITNESS_METRIC_HPP_

class quadratic_fitness_metric {
public:
    typedef double real_type;
    typedef real_type result_type;

    quadratic_fitness_metric( real_type s = 1. ) : m_scale(s) {}

    inline result_type operator()( real_type x ) {
        return operator()( x, m_scale );
    }

    inline result_type operator()( real_type x, real_type s ) {
        assert( s != 0 );

        result_type res = x / s;
        res *= res;

        res = 1.0 - res;

        if( res < 0.0 ) {
            return 0.0;
        }

        return res;
    }

    inline result_type operator()( const std::vector< real_type > & multi_variate ) {
        return operator()( multi_variate, m_scale );
    }

    inline result_type operator()( const std::vector< real_type > & multi_variate, real_type s ) {
        if( multi_variate.empty() ) return 1.0;

        return operator()( multi_variate.front(), s );
    }

protected:
    real_type m_scale;
};

#endif  // QUADRATIC_FITNESS_METRIC_HPP_
