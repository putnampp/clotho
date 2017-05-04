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
#ifndef IFITNESS_HPP_
#define IFITNESS_HPP_

#include <string>
#include <vector>
#include <ostream>

struct ifitness {
    typedef double result_type;

    /**
     * Reduce K-traits to a single 1D fitness value
     */
    virtual double operator()( const std::vector< double > & pheno ) = 0;
    virtual double operator()( double * first, double * last ) = 0;

    virtual double operator()( float * first, float * last ) = 0;
    virtual double operator()( const std::vector< float > & pheno ) = 0;

    virtual const std::string name() const = 0;

    virtual void log( std::ostream & out ) const = 0;

    virtual ~ifitness() {}
};

#endif  // IFITNESS_HPP_
