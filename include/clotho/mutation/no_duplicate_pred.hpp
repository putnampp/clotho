#ifndef CLOTHO_NO_DUPLICATE_HPP_
#define CLOTHO_NO_DUPLICATE_HPP_

namespace clotho {
namespace mutations {

template < class Set >
class no_duplicate_pred {
public:
    typedef Set                 set_type;
    typedef typename Set::element_type   element_type;

    no_duplicate_pred( set_type * elements ) :
        m_elements( elements )
    {}

    bool operator()( const element_type & elem ) {
        assert(false );
    }

    virtual ~no_duplicate_pred() {}
protected:
    set_type    * m_elements;
};

}   // namespace mutations {
}   // namespace clotho {

#endif  // CLOTHO_NO_DUPLICATE_HPP_
