#ifndef IFITNESS_GENERATOR_HPP_
#define IFITNESS_GENERATOR_HPP_

#include "clotho/genetics/ifitness.hpp"

#include <boost/property_tree/ptree.hpp>
#include <memory>

struct ifitness_generator {
    /**
     *  Creates a ifitness_generator of the derived type
     *  based upon the provided configuration parameters
     */
    virtual std::shared_ptr< ifitness_generator > create( boost::property_tree::ptree & config ) const = 0;

    /**
     * Generate the associated ifitness metric based upon
     * the provided population phenotype distribution
     */
    virtual std::shared_ptr< ifitness > generate( const std::vector< std::vector< double > > & ) = 0;
    virtual const std::string name() const = 0;

    virtual void log( std::ostream &  ) const = 0;

    virtual ~ifitness_generator() {}
};

#endif  // IFITNESS_GENERATOR_HPP_
