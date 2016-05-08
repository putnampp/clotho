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
#ifndef CLOTHO_RANDOM_SAMPLE_GENERATOR_HPP_
#define CLOTHO_RANDOM_SAMPLE_GENERATOR_HPP_

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class random_sample_generator {
public:
    typedef RNG                                 random_engine_type;
    typedef GeneticSpaceType                    genetic_space_type;
    typedef std::vector< size_t >                  sample_set_type;

    typedef typename sample_set_type::iterator sample_iterator;

    random_sample_generator( random_engine_type * rng, genetic_space_type * gs, size_t samp_size, bool with_replacement = false ) :
        m_rand(rng)
        , m_pop( gs )
    {
        generate( gs->sequence_count(), samp_size, with_replacement );
    }

    sample_iterator begin() {
        return m_sample.begin();
    }

    sample_iterator end() {
        return m_sample.end();
    }

    genetic_space_type * getPopulation() {
        return m_pop;
    }

    size_t getSampleAt( size_t idx ) {
        assert( 0 <= idx && idx < m_sample.size() );

        return m_sample[idx];
    }

    virtual ~random_sample_generator() {}
protected:
    void generate(size_t pop_size, size_t target_size, bool with_replacement ) {
        boost::random::uniform_int_distribution< size_t > dist( 0, pop_size - 1 );
        m_sample.reserve( target_size );

        if( with_replacement ) {
            while( m_sample.size() < target_size ) {
                size_t t = dist( *m_rand );
                m_sample.push_back( t );
            }
        } else if( target_size >= pop_size ) { // all samples should be unique
            size_t i = 0;
            while( i < target_size ) {
                m_sample.push_back( i++ );
            }
        } else {
            //std::vector< bool > indices( m_pop->sequence_count(), false );
            boost::dynamic_bitset< >    indices( m_pop->sequence_count(), false );

            while( m_sample.size() < target_size ) {
                size_t t = dist( *m_rand );
                while( indices[ t ] ) {
                    t = dist( *m_rand );
                }

                indices[ t ] = true;
                m_sample.push_back( t );
            }
        }
    }

    random_engine_type  * m_rand;
    genetic_space_type  * m_pop;
    sample_set_type     m_sample;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_RAtarget_sizeDOM_SAMPLE_GEtarget_sizeERATOR_HPP_
