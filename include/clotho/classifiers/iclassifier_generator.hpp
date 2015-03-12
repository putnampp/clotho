#ifndef ICLASSIFIER_GENERATOR_HPP_
#define ICLASSIFIER_GENERATOR_HPP_

#include "clotho/classifiers/iclassifier.hpp"

struct iclassifier_generator {

    virtual std::shared_ptr< iclassifier_generator > create( boost::property_tree::ptree & config ) = 0;

    virtual std::shared_ptr< iclassifier > generate( ) = 0;

    virtual const std::string name() const = 0;
    virtual void log( std::ostream & ) const = 0;

    virtual ~iclassifier_generator();
};

#endif  // ICLASSIFIER_GENERATOR_HPP_
