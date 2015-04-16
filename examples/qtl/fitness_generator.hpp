#ifndef FITNESS_GENERATOR_HPP_
#define FITNESS_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <iterator>
#include "clotho/genetics/normal_fitness_metric.hpp"

class normal_fitness_generator {
public:
    typedef normal_fitness_metric   result_type;

    normal_fitness_generator( boost::property_tree::ptree & config ) :
        m_mu(1.) {
        parseConfig( config );
    }

    template < class Iter >
    result_type operator()( Iter first, Iter last ) {
        // theoretical standard deviation:
        // sqrt( 2 * N * mu), where
        //  N - is the haploid sequence count
        //  mu - mutation rate per sequence
        // distance(,) - number of individuals in population (N_p) => N = 2 * (N_p); hence:  4.0 * N_p
        double n = 4.0 * (double)std::distance(first, last);

        n *= m_mu;
        n = sqrt(n);    // theoretical standard deviation

        return result_type( 0., n );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
        if( config.get_child_optional( oss.str() ) != boost::none ) {
            m_mu = config.get< double >( oss.str(), 1. );
        }
    }

    double m_mu;
};

#include "clotho/genetics/quadratic_fitness_metric.hpp"

extern const string FITNESS_BLOCK_K;
extern const string QUADRATIC_SCALE_K;

class quadratic_fitness_generator {
public:
    typedef quadratic_fitness_metric   result_type;

    quadratic_fitness_generator( boost::property_tree::ptree & config ) :
        m_scale(1.)
        , m_mu( 1. ) {
        parseConfig( config );
    }

    template < class Iter >
    result_type operator()( Iter first, Iter last ) {
        // theoretical standard deviation:
        // sqrt( 2 * N * mu), where
        //  N - is the haploid sequence count
        //  mu - mutation rate per sequence
        // distance(,) - number of individuals in population (N_p) => N = 2 * (N_p); hence:  4.0 * N_p
        double res = 4.0 * (double)std::distance( first, last );

        res *= m_mu;
        res = sqrt( res );  // theoretical standard deviation

        return result_type( m_scale * res );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << FITNESS_BLOCK_K << "." << QUADRATIC_SCALE_K;
        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_scale );
        } else {
            m_scale = config.get< double >( oss.str(), 1. );

            assert( m_scale > 0.0 );
        }

        oss.str("");
        oss.clear();
        oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
        if( config.get_child_optional( oss.str() ) != boost::none ) {
            m_mu = config.get< double >( oss.str(), 1. );
        }
    }

    double m_scale, m_mu;
};

#include "clotho/genetics/constant_fitness_metric.hpp"

extern const string CONSTANT_K;

class constant_fitness_generator {
public:
    typedef constant_fitness_metric   result_type;

    constant_fitness_generator( boost::property_tree::ptree & config ) :
        m_val( 1. ) {
        parseConfig( config );
    }

    template < class Iter >
    result_type operator()( Iter first, Iter last ) {
        return result_type( m_val );
    }
protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss << CONFIG_BLOCK_K << "." << FITNESS_BLOCK_K << "." << CONSTANT_K;
        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_val );
        } else {
            m_val = config.get< double >( oss.str(), 1. );
        }

    }

    double m_val;
};

#endif  // FITNESS_GENERATOR_HPP_
