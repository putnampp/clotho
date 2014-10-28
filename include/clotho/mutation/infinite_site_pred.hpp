#ifndef CLOTHO_INFINITE_SITE_PRED_HPP_
#define CLOTHO_INFINITE_SITE_PRED_HPP_

#include "clotho/mutation/no_duplicate_pred.hpp"

namespace clotho {
namespace mutations {

template < class Set >
struct infinite_site_pred {
    typedef no_duplicate_pred< Set > type;
};

}   // namespace mutations
}   // namespace clotho

#endif  // CLOTHO_INFINITE_SITE_PRED_HPP_
