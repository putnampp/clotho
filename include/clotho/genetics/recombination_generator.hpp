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
#ifndef RECOMBINATION_GENERATOR_HPP_
#define RECOMBINATION_GENERATOR_HPP_

#include "clotho/utility/random_generator.hpp"
#include "clotho/recombination/recombination_def.hpp"

#include <boost/random/bernoulli_distribution.hpp>

// Base Sequence Bias allows one to specify
// whether a parental sequence is more likely
// to be the base sequence to be passed along
// For example, the maternally inheritted sequence
// is more likely to be passed along than the
// paternal.
extern const string BASE_SEQUENCE_BIAS_K;

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Classifier, class T0, class T1 >
class random_generator< URNG, clotho::recombine::recombination< Sequence, Classifier, T0, T1 > > {
public:
    typedef URNG                                                                rng_type;
    typedef clotho::recombine::recombination< Sequence, Classifier, T0, T1 >    result_type;

    typedef clotho::utility::random_generator< URNG, Classifier >       classifier_generator_type;
    typedef typename classifier_generator_type::result_type             classifier_type;

    typedef double real_type;
    typedef boost::random::bernoulli_distribution< real_type >          dist_type;

    random_generator( rng_type & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_cgen( rng, config )
        , m_dist( 0.5 )
    {
        parseConfig( config );
    }

    random_generator( rng_type & rng, classifier_generator_type & cgen, real_type p = 0.5 ) :
        m_rng( &rng )
        , m_cgen( cgen )
        , m_dist( p ) 
    {    }

    result_type operator()() {
        classifier_type cfier = m_cgen();
        return result_type( cfier, m_dist( *m_rng ) );
    }

protected:

    void parseConfig( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( REC_BLOCK_K, lconfig );

        real_type p = lconfig.get< real_type >( BASE_SEQUENCE_BIAS_K, 0.5 );
        typename dist_type::param_type tmp(p);
        m_dist.param( tmp );

        lconfig.put( BASE_SEQUENCE_BIAS_K, p );
        config.put_child( REC_BLOCK_K, lconfig );
    }

    rng_type    * m_rng;
    classifier_generator_type m_cgen;
    dist_type   m_dist;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // RECOMBINATION_GENERATOR_HPP_
