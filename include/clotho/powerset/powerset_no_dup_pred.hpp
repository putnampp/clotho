#ifndef CLOTHO_POWERSET_NO_DUP_PRED_HPP_
#define CLOTHO_POWERSET_NO_DUP_PRED_HPP_

#include "clotho/powerset/powerset.hpp"
#include "clotho/mutation/no_duplicate_pred.hpp"

namespace clotho {
namespace mutations {

template < class Element, class Subset, class Block, class BlockMap, class ElementKeyer >
class no_duplicate_pred< clotho::powersets::powerset< Element, Subset, Block, BlockMap, ElementKeyer> > {
public:
    typedef clotho::powersets::powerset< Element, Subset, Block, BlockMap, ElementKeyer> set_type;
    typedef Element   element_type;

    no_duplicate_pred( set_type * elements ) :
        m_elements( elements ) {
    }

    bool operator()( const element_type & elem ) {
        std::pair< typename set_type::element_index_type, bool > res = m_elements->find_or_create( elem );
        return res.second;
    }

    virtual ~no_duplicate_pred() {}
protected:
    set_type    * m_elements;
};

}   // namespace mutations {
}   // namespace clotho {

#endif  // CLOTHO_POWERSET_NO_DUP_PRED_HPP_
