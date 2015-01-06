#ifndef REGION_CLASSIFIER_HPP_
#define REGION_CLASSIFIER_HPP_

#include <vector>
#include <algorithm>
#include <iterator>
#include <ostream>

namespace clotho {
namespace classifiers {

struct is_even_tag {};
struct is_odd_tag {};

/**
 *  A region based classifier.
 *
 *  Assumes that region_upper_bounds represent the upper bounds of each key region.
 *  For example, the regions of {0.25, 0.5, 0.75} would then be
 *  (-inf, 0.25) -> base
 *  [0.25, 0.5) -> alt
 *  [0.5, 0.75) -> base
 *  [0.75, inf) -> alt
 *
 *  The classifier expects a maps a specified element to an indexed region, and applies
 *  a Tag based ordering assignment.
 *
 *  The result is a bit mask which classifies the list elements at relative indices
 *  as being either in a "base" region (1), or "alt" region (0).
 */
template < class Element, class Result = bool, class Tag = is_even_tag >
class region_classifier {
public:
    typedef Element element_type;
    typedef Result  result_type;

    typedef std::vector< element_type > region_upper_bounds;
    typedef typename region_upper_bounds::const_iterator iterator;

    typedef region_upper_bounds     param_type;

    region_classifier() {}

    region_classifier( const param_type & bounds ) {
        reset_bounds( bounds );
    }

    region_classifier( const region_classifier< Element, Result, Tag > & rhs ) :
        m_ubounds( rhs.m_ubounds ) {
    }

    inline bool operator()( const element_type & elem ) const {
        if( m_ubounds.size() == 0 ) return 0;

        iterator it = std::upper_bound( m_ubounds.begin(), m_ubounds.end(), elem );

        return classify_region_index( it - m_ubounds.begin(), (Tag *)(NULL) );
    }

    template < class ElementIterator >
    bool operator()( ElementIterator it, size_t offset ) {
        return operator()( *(it + offset) );
    }

    void reset_bounds( const region_upper_bounds & bounds ) {
        m_ubounds.clear();

        if( std::is_sorted( bounds.begin(), bounds.end() ) ) {
            std::unique_copy( bounds.begin(), bounds.end(), std::back_inserter( m_ubounds ) );
        } else {
            region_upper_bounds tmp = bounds;

            std::sort( tmp.begin(), tmp.end() );
            std::unique_copy( tmp.begin(), tmp.end(), std::back_inserter( m_ubounds ) );
        }
    }

    size_t  region_count() const {
        return m_ubounds.size() + 1;
    }

    friend std::ostream & operator<<( std::ostream & os, const region_classifier & r ) {
        typename region_upper_bounds::const_iterator it = r.m_ubounds.begin();
        os << "{";
        if( it != r.m_ubounds.end() ) {
            os << *it;
            while( ++it != r.m_ubounds.end() ) {
                os << "," << *it;
            }
        }
        os << "}";
        return os;
    }

    virtual ~region_classifier() {
        m_ubounds.clear();
    }

protected:
    inline bool classify_region_index( unsigned int idx, is_even_tag * t ) const {
        return ((idx % 2) == 0);
    }

    inline bool classify_region_index( unsigned int idx, is_odd_tag * t ) const {
        return ((idx % 2) == 1);
    }

    region_upper_bounds m_ubounds;
};

}   // namespace classifier {
}   // namespace clotho {

#endif  // REGION_CLASSIFIER_HPP_
