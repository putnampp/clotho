#ifndef TEST_ELEMENT_H_
#define TEST_ELEMENT_H_

#include "clotho/powerset/block_map.hpp"
#include "clotho/powerset/element_key_of.hpp"
#include "clotho/powerset/normalized_key.hpp"

struct test_element {
    double k, v;

    test_element( double _k = 0., double _v = 0.) : k (_k), v(_v) {}
    test_element( const test_element & t ) : k(t.k), v(t.v) {}

    friend bool operator<( const test_element & lhs, const test_element & rhs );
    friend bool operator==( const test_element & lhs, const test_element & rhs ); 
    friend bool operator!=( const test_element & lhs, const test_element & rhs );
    friend std::ostream & operator<<( std::ostream & out, const test_element & rhs );
};

inline bool operator<( const test_element & lhs, const test_element & rhs ) {
    return lhs.k < rhs.k;
}

inline bool operator==( const test_element & lhs, const test_element & rhs ) {
    return ( lhs.k == rhs.k);
}

inline bool operator!=( const test_element & lhs, const test_element & rhs ) {
    return ( lhs.k != rhs.k);
}

inline std::ostream & operator<<( std::ostream & out, const test_element & rhs ) {
    out << "{" << rhs.k << ", " << rhs.v << "}";
    return out;
}

namespace clotho {
namespace powersets {

template <>
struct element_key_of< test_element > {
    typedef double key_type;

    inline key_type operator()( const test_element & t ) { return t.k; }

    static key_type get_key( const test_element & t ) { return t.k; }
};

template <>
struct normalized_key< test_element > : public key_range< > {
    inline double operator()( const test_element & elem ) {
        return elem.k;
    }

    static double get_key( const test_element & elem ) { return elem.k; }
};

}   // namespace powersets
}   // namespace clotho

#endif  // TEST_ELEMENT_H_
