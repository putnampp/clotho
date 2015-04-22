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
#ifndef IPOP_GROWTH_GENERATOR_HPP_
#define IPOP_GROWTH_GENERATOR_HPP_

#include "clotho/genetics/ipopulation_growth.hpp"

#include <boost/property_tree/ptree.hpp>
#include <memory>

struct ipopulation_growth_generator {
    /**
     *  Creates a ipopulation_growth_generator of the derived type
     *  based upon the provided configuration parameters
     */
    virtual std::shared_ptr< ipopulation_growth_generator > create( boost::property_tree::ptree & config ) const = 0;

    virtual std::shared_ptr< ipopulation_growth > generate() const = 0;

    virtual const std::string name() const = 0;

    virtual void log( std::ostream &  ) const = 0;

    virtual ~ipopulation_growth_generator() {}
};

#endif  // IPOP_GROWTH_GENERATOR_HPP_
