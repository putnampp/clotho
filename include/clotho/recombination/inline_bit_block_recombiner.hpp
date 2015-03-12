#ifndef INLINE_BIT_BLOCK_RECOMBINER_HPP_
#define INLINE_BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/recombination/bit_block_recombiner_def.hpp"
#include "clotho/recombination/inspect_methods.hpp"

#include "clotho/utility/bit_masks.hpp"

namespace clotho {
namespace recombine {
namespace walker {
namespace tag {

/**
 *  Inline classification walks each set bit in a block performing
 *  the classification step inline
 *
 *  This results in:
 *      - a single loop per 32-bit block (assuming 32 bit deBruijn hashing function)
 *      - uses a switch to jump between hashed keys
 *      - requires less bit shifting as mask offset bits can be inlined into case logic
 *  Empirically performs efficiently, though slightly less efficiently than earlier linear bit block iterator
 *  Still not as efficient as vector< index > approach.
 *
 *  Thought process:
 *      - a 'single' iteration per 32-bit block is more efficient than multiple
 *      - requires less memory as offsets are processed inline
 *
 *  Concerns:
 *      - Inline analysis of bit positions may result in less streamline execution (more context switching)
 */
struct inline_classify {};

}   // namespace tag
}   // namespace walker
}   // namespace recombine
}   // namespace clotho

namespace clotho {
namespace recombine {

#define CHECK_0()   if( m_cfier( *first ) ) { res |= (Block) OFFSET( 0 ); }
#define CHECK( x )  if( m_cfier( *(first + x) ) ) { res |= (Block) OFFSET( x ); }

template < class Classifier, class InspectMethodTag >
class bit_block_recombiner< Classifier, InspectMethodTag, clotho::recombine::walker::tag::inline_classify > {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( const Block b0, const Block b1, const ElementIterator first ) {
        Block b = InspectMethodTag::select(b0, b1);

        Block mask = walk(b, first);

        Block rec = ((b0 & mask) | (b1 & ~mask));

        return rec;
    }

protected:

    template < class ElementIterator >
    unsigned int walk( unsigned int lo, const ElementIterator first ) {
        typedef unsigned int Block;

        unsigned int res = 0;
        while( lo ) {
            unsigned int tmp = LEAST_SIG_BIT( lo );
            switch( DEBRUIJNBIT_HASH( tmp ) ) {
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
            lo ^= tmp;
        }
        return res;
    }

    template < class ElementIterator >
    unsigned long walk( unsigned long b, const ElementIterator first ) {
        unsigned long res = (unsigned long)0;
        unsigned int lo = (unsigned int) b;

        if( lo ) {
            res = (unsigned long) walk( lo, first );
        }

        lo = (unsigned int)( b >> 32 );

        if( lo ) {
            ElementIterator tmp = (first + 32);
            unsigned long hi_mask = (unsigned long)walk( lo, tmp );

            res |= (hi_mask << 32);
        }
        return res;
    }

    template < class ElementIterator >
    unsigned long unrolled_walk( unsigned long b, const ElementIterator first ) {
        typedef unsigned long Block;

        unsigned long res = (unsigned long)0;
        unsigned int lo = (unsigned int) b;

        while( lo ) {
            unsigned int tmp = LEAST_SIG_BIT( lo );
            switch( DEBRUIJNBIT_HASH( tmp ) ) {
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
            lo ^= tmp;
        }

        lo = (unsigned int) (b >> 32);
        while( lo ) {
            unsigned int tmp = LEAST_SIG_BIT( lo );
            switch( DEBRUIJNBIT_HASH( tmp ) ) {
            case 0:
                CHECK(32) break;    // 0
            case 1:
                CHECK(33) break;    // 1
            case 2:
                CHECK(60) break;    // 28
            case 3:
                CHECK(34) break;    // 2
            case 4:
                CHECK(61) break;    // 29
            case 5:
                CHECK(46) break;    // 14
            case 6:
                CHECK(56) break;    // 24
            case 7:
                CHECK(35) break;    // 3
            case 8:
                CHECK(62) break;    // 30
            case 9:
                CHECK(54) break;    // 22
            case 10:
                CHECK(52) break;    // 20
            case 11:
                CHECK(47) break;    // 15
            case 12:
                CHECK(57) break;    // 25
            case 13:
                CHECK(49) break;    // 17
            case 14:
                CHECK(36) break;    // 4
            case 15:
                CHECK(40) break;    // 8
            case 16:
                CHECK(63) break;    // 31
            case 17:
                CHECK(59) break;    // 27
            case 18:
                CHECK(45) break;    // 13
            case 19:
                CHECK(55) break;    // 23
            case 20:
                CHECK(53) break;    // 21
            case 21:
                CHECK(51) break;    // 19
            case 22:
                CHECK(48) break;    // 16
            case 23:
                CHECK(39) break;    // 7
            case 24:
                CHECK(58) break;    // 26
            case 25:
                CHECK(44) break;    // 12
            case 26:
                CHECK(50) break;    // 18
            case 27:
                CHECK(38) break;    //  6
            case 28:
                CHECK(43) break;    // 11
            case 29:
                CHECK(37) break;    //  5
            case 30:
                CHECK(42) break;    // 10
            case 31:
                CHECK(41) break;    //  9
            default:
                break;
            }
            lo ^= tmp;
        }
        return res;
    }

    classifier_type m_cfier;
};

#undef CHECK_0
#undef CHECK

}   // namespace recombine
}   // namespace clotho
#endif  // INLINE_BIT_BLOCK_RECOMBINER_HPP_
