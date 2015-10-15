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
#ifndef SEQUENCE_MUTATOR_GENERATOR_HPP_
#define SEQUENCE_MUTATOR_GENERATOR_HPP_

#include "clotho/utility/random_generator.hpp"
#include <boost/random/poisson_distribution.hpp>
#include "clotho/mutation/mutation_rate_parameter.hpp"

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Generator >
class random_generator< URNG, sequence_mutator< Sequence, Generator > > {
public:
    typedef URNG                                        rng_type;
    typedef sequence_mutator< Sequence, Generator >     result_type;

    typedef typename result_type::generator_result_type           allele_type;
    typedef typename allele_type::real_type             real_type;

    typedef mutation_rate_parameter< real_type >        mutation_type;
    typedef Generator                                   mutation_generator_type;

    typedef boost::random::poisson_distribution< unsigned int, real_type > dist_type;

//    random_generator( rng_type & rng,  mutation_generator_type & mgen, double rate ) : m_rng(&rng), m_mgen( mgen ), m_dist(rate) {}
    random_generator( rng_type & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_mut( config )
        , m_mgen( rng, config )
        , m_dist( mutation_type::DEFAULT_MUTATION_RATE ) {
        typename dist_type::param_type p( m_mut.m_mu );
        m_dist.param( p );
    }

    result_type operator()( unsigned int age = 0) {
        unsigned int nEvents = m_dist(*m_rng);
        return result_type( m_mgen, nEvents, age );
    }

protected:
    rng_type    * m_rng;
    mutation_type   m_mut;
    mutation_generator_type m_mgen;
    dist_type   m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // SEQUENCE_MUTATOR_GENERATOR_HPP_
