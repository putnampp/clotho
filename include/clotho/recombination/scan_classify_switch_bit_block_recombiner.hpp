#ifndef SCAN_CLASSIFY_SWITCH_BIT_BLOCK_RECOMBINER_HPP_
#define SCAN_CLASSIFY_SWITCH_BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/recombination/bit_block_recombiner_def.hpp"
#include "clotho/recombination/inspect_methods.hpp"

#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace recombine {
namespace walker {
namespace tag {

/**
 * Scan and Classify walks the set bits in a block and stores their
 * index (explodes the block into list of bit indices). The list of
 * indices is iterated and the elements at these offsets are classified.
 *
 * This results in:
 *      - a single bit walking loop per 32-bit block
 *      - a simple loop over offsets
 *      - switch over the bit indices
 *
 * Thought process:
 *      - Bit scanning could more streamline (less context switching)
 *      - Requires more memory (64 * 4 = 256 Bytes)
 *      - Simpler classification loop (indices will be ordered as a result of scan)
 *      - Switch logic allows mask bits sets to be defined inline rather than computed
 *
 *  Concern:
 *      - Doubles the number of iterative steps
 */
struct scan_and_classify_switch {};

}   // namespace tag
}   // namespace walker
}   // namespace recombine
}   // namespace clotho

namespace clotho {
namespace recombine {

#define CHECK_0()   if( m_cfier( *first ) ) { res |= (Block) OFFSET( 0 ); }
#define CHECK( x )  if( m_cfier( *(first + x) ) ) { res |= (Block) OFFSET( x ); }

template < class Classifier, class InspectMethodTag >
class bit_block_recombiner< Classifier, InspectMethodTag, clotho::recombine::walker::tag::scan_and_classify_switch > {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( const Block b0, const Block b1, const ElementIterator first ) {

        Block mask = walk( InspectMethodTag::select(b0, b1), first );

        Block rec = ((b0 & mask) | (b1 & ~mask));
        return rec;
    }

protected:
    classifier_type m_cfier;
    unsigned int indices[ 64 ];

    template < class ElementIterator >
    unsigned int walk( unsigned int b, const ElementIterator first ) {
        // Only compute the hash of each set bit
        // In effect, delay the offset lookup to be performed by the switching logic
        // Goal is to eliminate a 'double-lookup' (lookup from hash table, lookup based upon hash)
        // This assumes switch logic results in a 'jump table'
        //
        typedef unsigned int Block;
        unsigned int count = clotho::utility::hash_set_bits( b, indices );

        Block res = (Block)0;
        unsigned int * idx = indices;
        while( count-- ) {
            switch( *idx++ ) {
            case 0:
                CHECK_0() break;
            case 1:
                CHECK(1) break;
            case 2:
                CHECK(28) break;
            case 3:
                CHECK(2) break;
            case 4:
                CHECK(29) break;
            case 5:
                CHECK(14) break;
            case 6:
                CHECK(24) break;
            case 7:
                CHECK(3) break;
            case 8:
                CHECK(30) break;
            case 9:
                CHECK(22) break;
            case 10:
                CHECK(20) break;
            case 11:
                CHECK(15) break;
            case 12:
                CHECK(25) break;
            case 13:
                CHECK(17) break;
            case 14:
                CHECK(4) break;
            case 15:
                CHECK(8) break;
            case 16:
                CHECK(31) break;
            case 17:
                CHECK(27) break;
            case 18:
                CHECK(13) break;
            case 19:
                CHECK(23) break;
            case 20:
                CHECK(21) break;
            case 21:
                CHECK(19) break;
            case 22:
                CHECK(16) break;
            case 23:
                CHECK(7) break;
            case 24:
                CHECK(26) break;
            case 25:
                CHECK(12) break;
            case 26:
                CHECK(18) break;
            case 27:
                CHECK(6) break;
            case 28:
                CHECK(11) break;
            case 29:
                CHECK(5) break;
            case 30:
                CHECK(10) break;
            case 31:
                CHECK(9) break;
            default:
                break;
            }
        }
        return res;
    }

    template < class ElementIterator >
    unsigned long walk( unsigned long b, const ElementIterator first ) {
        // Only compute the hash of each set bit
        // In effect, delay the offset lookup to be performed by the switching logic
        // Goal is to eliminate a 'double-lookup' (lookup from hash table, lookup based upon hash)
        // This assumes switch logic results in a 'jump table'
        //

        unsigned int lo = (unsigned int)b;
        unsigned long mask = (unsigned long)walk( lo, first );

        lo = (unsigned int) ( b >> 32 );
        if( lo ) {
            ElementIterator tmp = (first + 32);
            unsigned long hi_mask = (unsigned long)walk(lo, tmp );
            mask |= (hi_mask << 32);
        }

        return mask;
    }

    template < class ElementIterator >
    unsigned long unrolled_walk( unsigned long b, const ElementIterator first ) {
        // Only compute the hash of each set bit
        // In effect, delay the offset lookup to be performed by the switching logic
        // Goal is to eliminate a 'double-lookup' (lookup from hash table, lookup based upon hash)
        // This assumes switch logic results in a 'jump table'
        //
        typedef unsigned long Block;
        unsigned int count = clotho::utility::hash_set_bits( b, indices );

        Block res = (Block)0;
        unsigned int * idx = indices;
        while( count-- ) {
            switch( *idx++ ) {
            case 0:
                CHECK_0() break;
            case 1:
                CHECK(1) break;
            case 2:
                CHECK(28) break;
            case 3:
                CHECK(2) break;
            case 4:
                CHECK(29) break;
            case 5:
                CHECK(14) break;
            case 6:
                CHECK(24) break;
            case 7:
                CHECK(3) break;
            case 8:
                CHECK(30) break;
            case 9:
                CHECK(22) break;
            case 10:
                CHECK(20) break;
            case 11:
                CHECK(15) break;
            case 12:
                CHECK(25) break;
            case 13:
                CHECK(17) break;
            case 14:
                CHECK(4) break;
            case 15:
                CHECK(8) break;
            case 16:
                CHECK(31) break;
            case 17:
                CHECK(27) break;
            case 18:
                CHECK(13) break;
            case 19:
                CHECK(23) break;
            case 20:
                CHECK(21) break;
            case 21:
                CHECK(19) break;
            case 22:
                CHECK(16) break;
            case 23:
                CHECK(7) break;
            case 24:
                CHECK(26) break;
            case 25:
                CHECK(12) break;
            case 26:
                CHECK(18) break;
            case 27:
                CHECK(6) break;
            case 28:
                CHECK(11) break;
            case 29:
                CHECK(5) break;
            case 30:
                CHECK(10) break;
            case 31:
                CHECK(9) break;
            case 32:
                CHECK(32) break;    // 0
            case 33:
                CHECK(33) break;    // 1
            case 34:
                CHECK(60) break;    // 28
            case 35:
                CHECK(34) break;    // 2
            case 36:
                CHECK(61) break;    // 29
            case 37:
                CHECK(46) break;    // 14
            case 38:
                CHECK(56) break;    // 24
            case 39:
                CHECK(35) break;    // 3
            case 40:
                CHECK(62) break;    // 30
            case 41:
                CHECK(54) break;    // 22
            case 42:
                CHECK(52) break;    // 20
            case 43:
                CHECK(47) break;    // 15
            case 44:
                CHECK(57) break;    // 25
            case 45:
                CHECK(49) break;    // 17
            case 46:
                CHECK(36) break;    // 4
            case 47:
                CHECK(40) break;    // 8
            case 48:
                CHECK(63) break;    // 31
            case 49:
                CHECK(59) break;    // 27
            case 50:
                CHECK(45) break;    // 13
            case 51:
                CHECK(55) break;    // 23
            case 52:
                CHECK(53) break;    // 21
            case 53:
                CHECK(51) break;    // 19
            case 54:
                CHECK(48) break;    // 16
            case 55:
                CHECK(39) break;    // 7
            case 56:
                CHECK(58) break;    // 26
            case 57:
                CHECK(44) break;    // 12
            case 58:
                CHECK(50) break;    // 18
            case 59:
                CHECK(38) break;    //  6
            case 60:
                CHECK(43) break;    // 11
            case 61:
                CHECK(37) break;    //  5
            case 62:
                CHECK(42) break;    // 10
            case 63:
                CHECK(41) break;    //  9
            default:
                break;
            }
        }
        return res;
    }
};

#undef CHECK_0
#undef CHECK

}   // namespace recombine
}   // namespace clotho
#endif  // SCAN_CLASSIFY_SWITCH_BIT_BLOCK_RECOMBINER_HPP_
