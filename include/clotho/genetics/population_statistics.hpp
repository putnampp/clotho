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
