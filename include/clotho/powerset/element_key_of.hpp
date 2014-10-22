#ifndef ELEMENT_KEY_OF_HPP_
#define ELEMENT_KEY_OF_HPP_


namespace clotho {
namespace powersets {

template < class E >
struct element_key_of {
    typedef E key_type;

    key_type operator()( const E & e ) {
        return e;
    }
};

}   // namespace powersets
}   // namespace cl

#endif  // ELEMENT_KEY_OF_HPP_
