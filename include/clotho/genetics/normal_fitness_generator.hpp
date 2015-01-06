#ifndef NORMAL_FITNESS_GENERATOR_HPP_
#define NORMAL_FITNESS_GENERATOR_HPP_

#include "clotho/genetics/ifitness_generator.hpp"
#include "clotho/genetics/normal_fitness_metric.hpp"

class normal_fitness_generator : public ifitness_generator {
public:
    typedef normal_fitness_metric   result_type;

    normal_fitness_generator();
    normal_fitness_generator( boost::property_tree::ptree & config );

    std::shared_ptr< ifitness_generator > create( boost::property_tree::ptree & config ) const;

    std::shared_ptr< ifitness > generate( const std::vector< std::vector< double > > & pop_traits );

    const std::string name() const;

    void log( std::ostream & out ) const;

    virtual ~normal_fitness_generator();

protected:

    void parseConfig( boost::property_tree::ptree & config );

    double m_mu;
};

#endif  // NORMAL_FITNESS_GENERATOR_HPP_
