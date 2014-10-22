#ifndef KEY_RANGE_HPP_
#define KEY_RANGE_HPP_

namespace clotho {
namespace powersets {

struct normal_01 {};

template < class R = normal_01 >
struct key_range {
    typedef R range_type;

    static constexpr double min = 0.0;
    static constexpr double max = 1.0;

    static constexpr double range = (max - min);
};

}   // namespace powersets
}   // namespace clotho

#endif  // KEY_RANGE_HPP_
