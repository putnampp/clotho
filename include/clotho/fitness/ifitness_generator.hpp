//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef IFITNESS_GENERATOR_HPP_
#define IFITNESS_GENERATOR_HPP_

#include "clotho/fitness/ifitness.hpp"

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
