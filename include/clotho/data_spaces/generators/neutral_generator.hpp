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
#ifndef CLOTHO_NEUTRAL_GENERATOR_HPP_
#define CLOTHO_NEUTRAL_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include <boost/random/bernoulli_distribution.hpp>

#include "clotho/data_spaces/generators/neutral_parameter.hpp"

namespace clotho {
namespace genetics {

template < class RNG >
class neutral_generator {
public:

    typedef neutral_generator< RNG >        self_type;

    typedef neutral_parameter< double >     parameter_type;
    typedef bool                            result_type;
    typedef RNG                             random_engine_type;

    neutral_generator( RNG * rng, const parameter_type & b ) :
        m_param( b )
        , m_rand( rng )
    {}

    neutral_generator( RNG * rng, boost::property_tree::ptree & config ) :
        m_param( config )
        , m_rand( rng )
    { }

    neutral_generator( const self_type & other ) :
        m_param( other.m_param )
        , m_rand( other.m_rand )
    {}

    result_type   operator()() {
        return m_param.m_dist( *m_rand );
    }

    virtual ~neutral_generator() {}
protected:
    parameter_type      m_param;
    random_engine_type  * m_rand;
};

class neutral_generator2 {
public:

    typedef neutral_parameter< double >     parameter_type;
    typedef bool                            result_type;

    neutral_generator2( boost::property_tree::ptree & config ) :
        m_param( config )
    {}

    neutral_generator2( const neutral_generator2 & other ) :
        m_param( other.m_param )
    {}

    template < class Engine >
    result_type operator()( Engine & eng ) {
        return m_param.m_dist( eng );
    }

    virtual ~neutral_generator2() {}
protected:
    parameter_type      m_param;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_NEUTRAL_GENERATOR_HPP_

