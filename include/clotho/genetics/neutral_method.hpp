#ifndef NEUTRAL_METHOD_HPP_
#define NEUTRAL_METHOD_HPP_

struct neutral_method {
    template < class E >
    static bool test( const E & elem ) {
        return true;
    }
};

#endif  // NEUTRAL_METHOD_HPP_
