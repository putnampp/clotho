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
#ifndef BASIC_ALLELE_GENERATOR_HPP_
#define BASIC_ALLELE_GENERATOR_HPP_

#include "basic_allele.h"
#include "clotho/utility/random_generator.hpp"
#include "clotho/utility/clotho_strings.hpp"

#include <boost/random/uniform_01.hpp>
#include <boost/random/bernoulli_distribution.hpp>

namespace clotho {
namespace utility {

template < class URNG >
class random_generator< URNG, basic_allele > {
public:
    typedef random_generator< URNG, basic_allele >  self_type;
    typedef basic_allele                            result_type;

    typedef typename result_type::real_type real_type;
    typedef typename result_type::key_type  key_type;

    typedef real_type selection_type;
    typedef real_type dominance_type;

    typedef boost::random::uniform_01< real_type >  dist_type;
    typedef boost::random::bernoulli_distribution< real_type > neutral_dist_type;

    random_generator( URNG & rng ) :
        m_rng( &rng )
        , m_neutral() {
    }

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_neutral() {
        parseConfig( config );
    }

    random_generator( const self_type & other ) :
        m_rng( other.m_rng )
        , m_neutral( other.m_neutral.param() ) {
    }

    result_type operator()( unsigned int age = 0 ) {
        basic_allele al;

        generate(al, age);

        return al;
    }

    virtual ~random_generator() {}
protected:

    void generate( basic_allele & a, unsigned int age ) {
        a.m_key = m_uniform( *m_rng );

        a.m_neutral = m_neutral( *m_rng );
        if( !a.m_neutral ) {
            a.m_select = m_uniform( *m_rng );
            a.m_dom = m_uniform( *m_rng );
        }
        a.m_age = age;
    }

    void parseConfig( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( ALLELE_BLOCK_K, lconfig );

        boost::property_tree::ptree nblock;
        nblock = lconfig.get_child( NEUTRAL_BLOCK_K, nblock );

        real_type p = nblock.get< real_type >( P_K, m_neutral.p() );

        nblock.put( P_K, p );

        lconfig.put_child( NEUTRAL_BLOCK_K, nblock );
        config.put_child( ALLELE_BLOCK_K, lconfig );

        typename neutral_dist_type::param_type tmp( p );
        m_neutral.param( tmp );
    }

    URNG *  m_rng;
    dist_type   m_uniform;
    neutral_dist_type m_neutral;
};

}   // namespace clotho {
}   // namespace utility {

#endif  // BASIC_ALLELE_GENERATOR_HPP_
