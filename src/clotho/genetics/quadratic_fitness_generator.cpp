#include "clotho/genetics/quadratic_fitness_generator.hpp"
#include "clotho/genetics/fitness_toolkit.hpp"

const std::string SCALE_K = "scale";
const std::string MU_K = "mu";

quadratic_fitness_generator::quadratic_fitness_generator() :
    m_scale( 1. )
    , m_mu( 1. ) {
    fitness_toolkit::getInstance()->register_tool( this );
}

quadratic_fitness_generator::quadratic_fitness_generator( boost::property_tree::ptree & config ) :
    m_scale(1.)
    , m_mu( 1. ) {
    parseConfig( config );
}

std::shared_ptr< ifitness_generator > quadratic_fitness_generator::create( boost::property_tree::ptree & config ) const {
    std::shared_ptr< ifitness_generator > t( new quadratic_fitness_generator( config ) );
    return t;
}

std::shared_ptr< ifitness > quadratic_fitness_generator::generate( const std::vector< std::vector< double > > & pop_traits ) {
    // theoretical standard deviation:
    // sqrt( 2 * N * mu), where
    //  N - is the haploid sequence count
    //  mu - mutation rate per sequence

    double res = 4.0 * (double)pop_traits.size();
    res *= m_mu;
    res = sqrt( res );  // theoretical standard deviation

    std::shared_ptr< ifitness > t( new result_type( m_scale * res ));
    return t;
}

const std::string quadratic_fitness_generator::name() const {
    return QUAD_NAME;
}

void quadratic_fitness_generator::parseConfig( boost::property_tree::ptree & config ) {
    std::ostringstream oss;
    oss /*<< QUAD_NAME << "."*/ << SCALE_K;
    if( config.get_child_optional( oss.str() ) == boost::none ) {
        config.put( oss.str(), m_scale );
    } else {
        m_scale = config.get< double >( oss.str(), m_scale );

        assert( m_scale > 0.0 );
    }

    oss.str("");
    oss.clear();
    oss /*<< QUAD_NAME << "."*/ << MU_K;
    if( config.get_child_optional( oss.str() ) == boost::none ) {
        config.put( oss.str(), m_mu );
    } else {
        m_mu = config.get< double >( oss.str(), m_mu );
    }
}

void quadratic_fitness_generator::log( std::ostream & out ) const {
    out << "{" << QUAD_NAME << "_gen, " << m_scale << ", " << m_mu << "}\n";
}

quadratic_fitness_generator::~quadratic_fitness_generator() {}

static const quadratic_fitness_generator qfg_reg;
