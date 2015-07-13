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
#ifndef ICONFIGURABLE_HPP_
#define ICONFIGURABLE_HPP_

#include <string>
#include <ostream>
#include <boost/program_options.hpp>

#include "commandline_actions.hpp"

namespace po=boost::program_options;

struct iconfigurable {

    virtual std::string name() const = 0;
    virtual void getOptions( po::options_description & od ) = 0;
    virtual COMMANDLINE_ACTION validate( const po::variables_map & vm ) = 0;

    virtual void printMessage( std::ostream & out ) {}
    virtual void printWarning( std::ostream & out ) {}

    virtual ~iconfigurable() {}
};

#endif  // ICONFIGURABLE_HPP_
