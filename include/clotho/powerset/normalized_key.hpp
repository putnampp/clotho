#ifndef NORMALIZED_KEY_HPP_
#define NORMALIZED_KEY_HPP_

#include "clotho/powerset/key_range.hpp"

namespace clotho {
namespace powersets {

template < class Element >
struct normalized_key : public key_range < Element > {
    inline double operator()( const Element & elem ) {
        return 0.0;
    }
};

}   // namespace powersets
}   // namespace clotho

#endif  // NORMALIZED_KEY_HPP_
