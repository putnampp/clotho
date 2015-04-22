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
#ifndef POPULATION_STATISTICS_HPP_
#define POPULATION_STATISTICS_HPP_

#include <unordered_map>
#include <vector>

template < class Element >
class population_statistics {
public:
    typedef Element                                     element_type;
    typedef std::unordered_map< element_type, size_t >  element_map_type;

    population_statistics();

    void add( const element_type & e );

    virtual ~population_statistics();

protected:
    element_map_type    m_elements;
};

template < class Element >
class population_statistics< std::shared_ptr< Element > > {
public:
    typedef std::shared_ptr< Element > element_type;
    typedef std::unordered_map< element_type, size_t >  element_map_type;

    typedef std::pair< Element, unsigned int >          element_frequency;
    typedef std::vector< element_frequency >            element_set;

    typedef element_map_type::iterator  element_iterator;

    void add( element_type elem ) {
        element_iterator it = m_element_lookup.find( elem );

        if( it == m_elements.end() ) {
            m_element_lookup.insert( std::make_pair( elem, m_elem_freq.size() ) );
            m_elem_freq.push_back( std::make_pair( elem, 1 ) );
        } else {
            ++(m_elem_freq[it->second].second);
        }
    }

    virtual ~population_statistics() {
        m_elem_freq.clear();
        m_elements.clear();
    }
protected:
    element_map_type        m_element_lookup;
    element_frequency_set   m_elem_freq;
};

#endif  // POPULATION_STATISTICS_HPP_
