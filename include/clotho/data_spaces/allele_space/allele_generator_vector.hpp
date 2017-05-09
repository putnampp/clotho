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
#ifndef CLOTHO_ALLELE_GENERATOR_VECTOR_HPP_
#define CLOTHO_ALLELE_GENERATOR_VECTOR_HPP_

#include "clotho/data_spaces/allele_space/allele_generator.hpp"
#include "clotho/data_spaces/allele_space/allele_space_vector.hpp"

#include "clotho/data_spaces/generators/position_generator.hpp"
//#include "clotho/data_spaces/generators/neutral_generator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PositionType, class SizeType >
class AlleleGenerator< RNG, AlleleSpace< PositionType, SizeType > > {
public:
    typedef AlleleGenerator< RNG, AlleleSpace< PositionType, SizeType > > self_type;

    typedef RNG                                     random_engine_type;

    typedef AlleleSpace< PositionType, SizeType >   allele_type;
    typedef typename allele_type::position_type     position_type;

    typedef position_generator2< position_type > position_generator_type;
//    typedef neutral_generator2                   neutrality_generator_type;

    AlleleGenerator( boost::property_tree::ptree & config ) :
        m_rng( NULL )
//        , m_neut_gen( config )
    {}

    AlleleGenerator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rng( rng )
//        , m_neut_gen( config )
    {}

    AlleleGenerator( random_engine_type * rng, const self_type & other ) :
        m_rng( rng )
//        , m_neut_gen( other.m_neut_gen )
    {}

    AlleleGenerator( const self_type & other ) :
        m_rng( other.m_rng )
//        , m_neut_gen( other.m_neut_gen )
    {}

    template < class FreeSpaceType >
    void operator()( allele_type & all, FreeSpaceType & free_space, unsigned int N, unsigned int age ) {
        generate( all, free_space, N, age );
    }

    template < class FreeSpaceType >
    void generate( allele_type & all, FreeSpaceType & free_space, unsigned int N, unsigned int age ) {
        assert( N <= free_space.free_size());

        typename FreeSpaceType::base_type::iterator it = free_space.free_begin();
        
        while( N-- ) {
            typename FreeSpaceType::size_type offset = *it++;
            operator()( all, offset, age );
        }
    }

    void operator()( allele_type & all, unsigned int idx, unsigned int age ) {
        position_type p = m_pos_gen( *m_rng );

//        bool neu = m_neut_gen( *m_rng );
//        all.setAllele( idx, p, neu, age );
//
        all.setAllele( idx, p, age );
    }

    virtual ~AlleleGenerator() {}

protected:
    random_engine_type          * m_rng;
    position_generator_type     m_pos_gen;
//    neutrality_generator_type   m_neut_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ALLELE_GENERATOR_VECTOR_HPP_
