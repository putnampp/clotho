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

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Generator >
class random_generator< URNG, sequence_mutator< Sequence, Generator > > {
public:
    typedef URNG                                        rng_type;
    typedef sequence_mutator< Sequence, Generator >     result_type;

    typedef Generator                                   mutation_generator_type;

    typedef boost::random::poisson_distribution< unsigned int, double > dist_type;

//    random_generator( rng_type & rng,  mutation_generator_type & mgen, double rate ) : m_rng(&rng), m_mgen( mgen ), m_dist(rate) {}
    random_generator( rng_type & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_mgen( rng, config )
        , m_dist( DEFAULT_MUTATION_RATE ) {
        parseConfig( config );
    }

    result_type operator()( unsigned int age = 0) {
        unsigned int nEvents = m_dist(*m_rng);
        return result_type( m_mgen, nEvents, age );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss /*<< CONFIG_BLOCK_K << "."*/ << MUT_BLOCK_K << "." << RATE_PER_REGION_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_dist.mean() );
        } else {
            double mean = config.get< double >( oss.str(), DEFAULT_MUTATION_RATE );
            typename dist_type::param_type p( mean );
            m_dist.param( p );
        }
    }

    rng_type    * m_rng;
    mutation_generator_type m_mgen;
    dist_type   m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // SEQUENCE_MUTATOR_GENERATOR_HPP_
