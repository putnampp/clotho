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
#ifndef TRAIT_WEIGHT_GENERATOR_HPP_
#define TRAIT_WEIGHT_GENERATOR_HPP_

#include "trait_weight.hpp"

#include <boost/random/normal_distribution.hpp>

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, trait_weight< double > > {
public:
    typedef double real_type;
    typedef real_type   result_type;
    typedef random_generator< URNG, trait_weight< real_type > > self_type;

    typedef boost::random::normal_distribution< real_type >    dist_type;

    class param_type {
    public:
        real_type          _mean, _sigma;

        param_type( real_type m = 0., real_type s = 1. ) :
            _mean(m), _sigma(s) {
        }

        double mean()   const {
            return _mean;
        }
        double sigma()  const {
            return _sigma;
        }
    };

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_dist( ) {
        parseConfig( config );
    }

    random_generator( URNG & rng, real_type mean = 0.0, real_type sigma = 1.0 ) :
        m_rng( &rng )
        , m_dist( mean, sigma ) {
    }

    random_generator( const self_type & other ) :
        m_rng( other.m_rng )
        , m_dist( other.m_dist.param() ) {
    }

    result_type operator()() {
        return m_dist(*m_rng);
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        std::ostringstream oss;
        oss /*<< CONFIG_BLOCK_K << "."*/ << TRAIT_BLOCK_K << "." << MEAN_K;

        real_type m = m_dist.mean(), s = m_dist.sigma();
        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), m_dist.mean() );
        } else {
            m = config.get< real_type >( oss.str(), 0. );
        }

        oss.str("");
        oss.clear();

        oss /*<< CONFIG_BLOCK_K << "."*/ << TRAIT_BLOCK_K << "." << SIGMA_K;

        if( config.get_child_optional( oss.str() ) == boost::none ) {
            config.put( oss.str(), s );
        } else {
            s = config.get< real_type >( oss.str(), 1. );
        }

        typename dist_type::param_type p( m, s );
        m_dist.param( p );
    }

    URNG            * m_rng;
    dist_type       m_dist;
};

}   // namespace utility
}   // namespace clotho

#endif  // TRAIT_WEIGHT_GENERATOR_HPP_
