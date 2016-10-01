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
#ifndef CLOTHO_GENERAL_FITNESS_SELECTION_GENERATOR_HPP_
#define CLOTHO_GENERAL_FITNESS_SELECTION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "clotho/data_spaces/selection/fitness_selection_generator.hpp"
#include "clotho/data_spaces/fitness/general_fitness.hpp"

namespace clotho {
namespace genetics {

template < class RNG >
class SelectionGenerator< RNG, fitness_selection< GeneralFitness > > {
public:
    typedef RNG             random_engine_type;
    typedef GeneralFitness  fitness_space_type;

    typedef typename fitness_space_type::individual_id_type                     individual_id_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_vector;

    typedef typename mate_pair_vector::iterator                                 iterator;
    typedef typename mate_pair_vector::const_iterator                           const_iterator;

    SelectionGenerator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
    {}

    void update( const fitness_space_type & fit, unsigned int N ) {
        typedef typename fitness_space_type::fitness_type   fitness_score_type;
        typedef typename fitness_space_type::const_iterator iterator;

        bool constant_fitness = true;
        unsigned int M = fit.individual_count();

        iterator it = fit.begin(), end = fit.end();

        if( it != end ) {
            fitness_score_type prev = *it++;
            while( constant_fitness && it != end ) {
                fitness_score_type cur = *it++;
                constant_fitness = (prev == cur);
            }
        }

        if( constant_fitness ) {
            boost::random::uniform_int_distribution< individual_id_type > uni( 0, ((M == 0) ? 0 : (M - 1)) );
            generate( uni, N );
        } else {
            boost::random::discrete_distribution< individual_id_type, fitness_score_type > disc( fit.begin(), fit.end() );
            generate( disc, N );
        }
    }

    mate_pair_vector & getMatePairs() const {
        return m_pairs;
    }
    
    unsigned int individual_count() const {
        return m_pairs.size();
    }

    iterator    begin() {
        return m_pairs.begin();
    }

    iterator    end() {
        return m_pairs.end();
    }

    const_iterator    begin() const {
        return m_pairs.begin();
    }

    const_iterator    end() const {
        return m_pairs.end();
    }

    unsigned int size() const {
        return individual_count();
    }

    virtual ~SelectionGenerator() {}

protected:

    template < class DistributionType >
    void generate( DistributionType & dist, size_t count ) {
        m_pairs.clear();

        while( m_pairs.size() < count ) {
            individual_id_type  id0 = dist( *m_rand );
            individual_id_type  id1 = dist( *m_rand );

            m_pairs.push_back( std::make_pair( id0, id1 ) );
        }
    }

    random_engine_type  * m_rand;
    mate_pair_vector    m_pairs;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_GENERAL_FITNESS_SELECTION_GENERATOR_HPP_

