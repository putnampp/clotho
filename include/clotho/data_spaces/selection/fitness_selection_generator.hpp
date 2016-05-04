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
#ifndef CLOTHO_FITNESS_SELECTION_GENERATOR_HPP_
#define CLOTHO_FITNESS_SELECTION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include <boost/random/discrete_distribution.hpp>

#include <vector>

#include "clotho/data_spaces/selection/selection_generator.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
struct fitness_selection {};

template < class RNG, class GeneticSpaceType >
class SelectionGenerator< RNG, fitness_selection< GeneticSpaceType > > {
public:
    typedef RNG                                                                 random_engine_type;
    typedef GeneticSpaceType                                                    genetic_space_type;
    typedef typename genetic_space_type::individual_id_type                     individual_id_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_vector;
    typedef typename mate_pair_vector::iterator                                 parent_iterator;
    typedef typename mate_pair_vector::const_iterator                           const_parent_iterator;

    typedef typename genetic_space_type::fitness_score_type                     fitness_score_type;
    typedef typename genetic_space_type::fitness_scores                         fitness_scores;

    typedef boost::random::discrete_distribution< individual_id_type, fitness_score_type >  distribution_type;

    SelectionGenerator( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_rand( rng )
    { }

    void update( genetic_space_type * parents, unsigned int count ) {
        typename distribution_type::param_type  param_type;
        m_dist.param( param_type( parents->fitness_begin(), parents->fitness_end() ) );

        m_pairs.clear();

        while( count-- ) {
            individual_id_type  id0 = m_dist( *m_rand );
            individual_id_type  id1 = m_dist( *m_rand );

            m_pairs.push_back( std::make_pair( id0, id1 ) );
        }   
    }

    parent_iterator begin() {
        return m_pairs.begin();
    }

    parent_iterator end() {
        return m_pairs.end();
    }

    const_parent_iterator begin() const {
        return m_pairs.begin();
    }

    const_parent_iterator end() const {
        return m_pairs.end();
    }

    mate_pair_vector & getMatePairs() {
        return m_pairs;
    }

    size_t size() const {
        return m_pairs.size();
    }

    virtual ~SelectionGenerator() {}

protected:
    random_engine_type  * m_rand;
    mate_pair_vector    m_pairs;

    distribution_type   m_dist;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_FITNESS_SELECTION_GENERATOR_HPP_
