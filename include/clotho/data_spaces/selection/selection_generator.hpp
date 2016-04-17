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
#ifndef CLOTHO_SELECTION_GENERATOR_HPP_
#define CLOTHO_SELECTION_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include <vector>

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
class SelectionGenerator {
public:
    typedef GeneticSpaceType    genetic_space_type;

    typedef std::vector< std::pair< size_t, size_t > > mate_pair_vector;
    typedef typename mate_pair_vector::iterator     parent_iterator;
    typedef typename mate_pair_vector::const_iterator     const_parent_iterator;

    SelectionGenerator( boost::property_tree::ptree & config ){}

    void update( genetic_space_type * parents, fitness_type *, unsigned int

    parent_iterator begin() {
        return m_pairs.begin();
    }

    parent_iterator end() {
        return m_pairs.end()
    }

    const_parent_iterator begin() const {
        return m_pairs.begin();
    }

    const_parent_iterator end() const {
        return m_pairs.end()
    }

protected:
    mate_pair_vector    m_pairs;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_SELECTION_GENERATOR_HPP_
