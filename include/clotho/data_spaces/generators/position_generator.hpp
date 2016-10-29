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
#ifndef CLOTHO_POSITION_GENERATOR_HPP_
#define CLOTHO_POSITION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/data_spaces/generators/position_distribution_helper.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PositionType >
class position_generator {
public:

    typedef position_generator< RNG, PositionType > self_type;

    typedef RNG                 random_engine_type;
    typedef PositionType        position_type;

    typedef typename position_distribution_helper< PositionType >::type distribution_type;

    position_generator( random_engine_type * rng ) :
        m_rand( rng )
    {}

    position_generator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
    {  }

    position_generator( const self_type & other ) :
        m_rand( other.m_rand )
    {}


    position_type   operator()() {
        return m_dist( *m_rand );
    }

    virtual ~position_generator() {}
protected:
    random_engine_type  * m_rand;
    distribution_type   m_dist;
};

template < class PositionType >
class position_generator2 {
public:
    typedef position_generator2< PositionType > self_type;
    typedef PositionType        position_type;
    typedef typename position_distribution_helper< PositionType >::type distribution_type;

    position_generator2() {}

    template < class Engine >
    position_type operator()( Engine & eng ) {
        return m_dist( eng );
    }

    virtual ~position_generator2() {}

protected:
    distribution_type   m_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POSITION_GENERATOR_HPP_
