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
#ifndef IPOP_GROWTH_HPP_
#define IPOP_GROWTH_HPP_

#include <string>
#include <vector>
#include <ostream>

struct ipopulation_growth {
/**
 *  Generates a population size based upon a population size (psize) or a generation (gen)
 *
 *  Both psize and gen are context specific.
 */
    virtual unsigned int operator()( unsigned int psize, unsigned int gen ) = 0;
    virtual const std::string name() const = 0;

    virtual void log( std::ostream & out ) const = 0;

    virtual ~ipopulation_growth() {}
};

#endif  // IPOP_GROWTH_HPP_

