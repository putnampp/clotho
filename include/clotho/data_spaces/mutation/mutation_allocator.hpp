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
#ifndef CLOTHO_MUTATION_ALLOCATOR_HPP_
#define CLOTHO_MUTATION_ALLOCATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/random/poisson_distribution.hpp>

#include "clotho/mutation/mutation_rate_parameter.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class IntType = size_t >
class mutation_allocator {
public:
    typedef RNG                                     random_engine_type;
    typedef IntType                                 int_type;
    typedef double                                  real_type;
    typedef mutation_rate_parameter< real_type >    rate_type;

    typedef boost::random::poisson_distribution< int_type, real_type > distribution_type;
    

    mutation_allocator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
        , m_mut_rate( config )
    {}

    int_type allocate( size_t N ) {
        typedef typename distribution_type::param_type param_type;

        m_dist.param( param_type( m_mut_rate.m_mu * N ) );

        return m_dist( *m_rand );
    }

    virtual ~mutation_allocator() {}

protected:
    random_engine_type  * m_rand;
    rate_type           m_mut_rate;
    distribution_type   m_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_MUTATION_ALLOCATOR_HPP_
