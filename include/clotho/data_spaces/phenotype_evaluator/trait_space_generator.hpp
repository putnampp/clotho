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
#ifndef CLOTHO_TRAIT_SPACE_GENERATOR_HPP_
#define CLOTHO_TRAIT_SPACE_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/generators/weight_parameter.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class TraitSpaceType >
class TraitSpaceGenerator {
public:
    typedef TraitSpaceGenerator< RNG, TraitSpaceType > self_type;

    typedef RNG                                     random_engine_type;
    typedef TraitSpaceType                          trait_space_type;
    typedef typename trait_space_type::weight_type  weight_type;

    typedef weight_parameter< weight_type >         parameter_type;

    TraitSpaceGenerator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
        , m_weights( config )
    {}

    TraitSpaceGenerator( random_engine_type * rng, const parameter_type & param ) :
        m_rng( rng )
        , m_weights( param )
    {}

    TraitSpaceGenerator( const self_type & other ) :
        m_rng( other.m_rng )
        , m_weights( other.m_weights )
    {}

    template < class FreeSpaceType >
    void operator()( trait_space_type & traits, FreeSpaceType & fs, unsigned int N ) {
        generate( traits, fs, N );
    }

    template < class FreeSpaceType >
    void generate( trait_space_type & traits, FreeSpaceType & fs, unsigned int N ) {
        assert( N <= fs.free_size() );

        typename FreeSpaceType::base_type::iterator it = fs.free_begin();

        while( N-- ) {
            typename FreeSpaceType::size_type offset = *it++;
            operator()( traits, offset );
        }
    }

    void operator()( trait_space_type & traits, unsigned int offset ) {
        unsigned int i = 0;
        while( i < traits.trait_count() ) {
            traits.setWeight( m_weights.m_dist( *m_rng ), offset, i );
            ++i;
        }
    }

    virtual ~TraitSpaceGenerator() {}

protected:
    random_engine_type  * m_rng;
    parameter_type      m_weights;
};

template < class TraitSpaceType >
class TraitSpaceGenerator2 {
public:

    typedef TraitSpaceGenerator2< TraitSpaceType >  self_type;
    typedef TraitSpaceType                          trait_space_type;
    typedef typename trait_space_type::weight_type  weight_type;

    typedef weight_parameter< weight_type >         parameter_type;

    TraitSpaceGenerator2( boost::property_tree::ptree & config ) :
        m_weights( config )
    {}

    TraitSpaceGenerator2( const parameter_type & param ) :
        m_weights( param )
    {}
    
    TraitSpaceGenerator2( const self_type & other ) :
        m_weights( other.m_weights )
    {}
    
    template < class Engine >
    void operator()( Engine & eng, trait_space_type & traits, unsigned int offset ) {
        unsigned int i = 0;
        while( i < traits.trait_count() ) {
            weight_type w = m_weights.m_dist( eng );
            traits.setWeight( w, offset, i );
            ++i;
        }
    }

    virtual ~TraitSpaceGenerator2() {}
protected:

    parameter_type m_weights;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_TRAIT_SPACE_GENERATOR_HPP_
