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
