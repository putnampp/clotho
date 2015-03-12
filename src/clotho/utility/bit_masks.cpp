#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace utility {

unsigned int scan_set_bits( unsigned int b, unsigned int * indices ) {
    unsigned int count = 0;
    while( b ) {
        unsigned int tmp = LEAST_SIG_BIT( b );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP( tmp );
        b ^= tmp;
    }
    return count;
}

unsigned int scan_set_bits( unsigned long b, unsigned int * indices ) {
    unsigned int count = 0;
    unsigned int lo = (unsigned int) b;
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP( tmp );
        lo ^= tmp;
    }

    lo = (unsigned int)( b >> 32 );
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH_LOOKUP_HI( tmp );
        lo ^= tmp;
    }
    return count;
}

unsigned int hash_set_bits( unsigned int b, unsigned int * indices ) {
    unsigned int count = 0;
    while( b ) {
        unsigned int tmp = LEAST_SIG_BIT( b );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp );
        b ^= tmp;
    }
    return count;
}

unsigned int hash_set_bits( unsigned long b, unsigned int * indices ) {
    unsigned int count = 0;
    unsigned int lo = (unsigned int) b;
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp );
        lo ^= tmp;
    }

    lo = (unsigned int)( b >> 32 );
    while( lo ) {
        unsigned int tmp = LEAST_SIG_BIT( lo );
        indices[ count++ ] = DEBRUIJNBIT_HASH( tmp ) + 32;
        lo ^= tmp;
    }
    return count;
}

}   // namespace utility {
}   // namespace clotho {
