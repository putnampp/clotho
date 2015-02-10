#ifndef RANDOM_FITNESS_METRIC_HPP_
#define RANDOM_FITNESS_METRIC_HPP_

#include "clotho/genetics/ifitness.hpp"

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
