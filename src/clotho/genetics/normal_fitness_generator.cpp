#include "clotho/genetics/normal_fitness_generator.hpp"
#include "clotho/genetics/fitness_toolkit.hpp"

const std::string MU_K = "mu";

normal_fitness_generator::normal_fitness_generator() : m_mu(1.) {
    fitness_toolkit::getInstance()->register_tool( this );
}

normal_fitness_generator::normal_fitness_generator( boost::property_tree::ptree & config ) :
    m_mu(1.) {
    parseConfig( config );
}

std::shared_ptr< ifitness_generator > normal_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new normal_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > normal_fitness_generator::generate( const std::vector< std::vector< double > > & pop_traits ) {
    // theoretical standard deviation:
    // sqrt( 2 * P * N * mu), where
    //  N - is the individual count
    //  P - is the ploidy of an individual = 2.0
    //  mu - mutation rate per sequence
    //
    //  NP = haploid sequence count of the population

    double n = 4.0 * (double)pop_traits.size();

    n *= m_mu;
    n = sqrt(n);    // theoretical standard deviation

    std::shared_ptr< ifitness > r( new result_type( 0., n ) );
    return r;
}

void normal_fitness_generator::parseConfig( boost::property_tree::ptree & config ) {
    std::ostringstream oss;
    oss /*<< NORM_NAME << "."*/ << MU_K;
    if( config.get_child_optional( oss.str() ) == boost::none ) {
        config.put( oss.str(), m_mu );
    } else {
        m_mu = config.get< double >( oss.str(), 1. );
    }
}

const std::string normal_fitness_generator::name() const {
    return NORM_NAME;
}

void normal_fitness_generator::log( std::ostream & out ) const {
    out << "{" << NORM_NAME << "_gen"
        << ", " << m_mu
        << "}\n";
}

normal_fitness_generator::~normal_fitness_generator() {}

static const normal_fitness_generator nfg_reg;
