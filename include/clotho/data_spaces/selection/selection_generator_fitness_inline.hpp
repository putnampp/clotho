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
#ifndef CLOTHO_SELECTION_GENERATOR_FITNES_INLINE_HPP_
#define CLOTHO_SELECTION_GENERATOR_FITNES_INLINE_HPP_

#include <boost/property_tree/ptree.hpp>

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "clotho/data_spaces/selection/selection_generator.hpp"
#include "clotho/data_spaces/fitness/general_fitness.hpp"

#include "clotho/data_spaces/selection/selection_details.hpp"

namespace clotho {
namespace genetics {

template < class RNG >
class UniformPairGenerator : public selection_details< RNG, unsigned int > {
public:
    typedef selection_details< RNG, unsigned int >    base_type;
    typedef typename base_type::individual_id_type individual_id_type;
    typedef typename base_type::parent_pair parent_pair;

    UniformPairGenerator( RNG * rng, unsigned int N ) :
        base_type(rng)
        , m_dist( 0, N)
    {  }

    parent_pair operator()() {
        return std::make_pair( nextID(), nextID() );
    }

    individual_id_type  nextID() {
        return m_dist( *(this->m_rand) );
    }

    virtual ~UniformPairGenerator() {}

protected:
    boost::random::uniform_int_distribution< individual_id_type > m_dist;
};

template < class RNG, class ScoreType >
class DiscretePairGenerator : public selection_details< RNG, unsigned int > {
public:

    typedef selection_details< RNG, unsigned int > base_type;
    typedef typename base_type::individual_id_type  individual_id_type;
    typedef typename base_type::parent_pair         parent_pair;
    typedef ScoreType                               score_type;

    template< class ScoreIterator >
    DiscretePairGenerator( RNG * rng, ScoreIterator first, ScoreIterator last ) :
        base_type( rng )
        , m_dist( first, last )
    {}

    parent_pair operator()() {
        return std::make_pair( nextID(), nextID() );
    }

    individual_id_type nextID() {
        return m_dist( *(this->m_rand) );
    }

    virtual ~DiscretePairGenerator() {}

protected:
    boost::random::discrete_distribution< individual_id_type, score_type > m_dist;   
};

template < class RNG >
class ConstantPairGenerator : public selection_details< RNG, unsigned int > {
public:

    typedef selection_details< RNG, unsigned int >  base_type;
    typedef typename base_type::individual_id_type  individual_id_type;
    typedef typename base_type::parent_pair         parent_pair;

    ConstantPairGenerator( RNG * rng, unsigned int id = 0 ) :
        base_type( rng )
        , m_id(id)
    {}

    parent_pair operator()() {
        return std::make_pair( m_id, m_id );
    }


    virtual ~ConstantPairGenerator() {}

protected:
    unsigned int m_id;
};

}   // namepsace_genetics
}   // namespace clotho

#endif  // CLOTHO_SELECTION_GENERATOR_FITNES_INLINE_HPP_

