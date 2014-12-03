#ifndef QTL_WEIGHT_ITERATOR_HPP_
#define QTL_WEIGHT_ITERATOR_HPP_

#include "clotho/utility/iterator_helper.hpp"
#include "qtl_allele.h"

namespace clotho {
namespace utility {

template < >
struct iterator_helper< qtl_allele > {
    typedef qtl_allele::weight_iterator type;

    static type make_first( qtl_allele & all ) {
        return all.begin();
    }

    static type make_last( qtl_allele & all ) {
        return all.end();
    }
};

}   // namespace utility {
}   // namespace clotho {
#endif  // QTL_WEIGHT_ITERATOR_HPP_
