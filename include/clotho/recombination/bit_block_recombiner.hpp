#ifndef BIT_BLOCK_RECOMBINER_HPP_
#define BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/recombination/iterate_bit_block_recombiner.hpp"
#include "clotho/recombination/scan_classify_bit_block_recombiner.hpp"
#include "clotho/recombination/scan_classify_switch_bit_block_recombiner.hpp"
#include "clotho/recombination/inline_bit_block_recombiner.hpp"

//#if defined( USE_SCAN_AND_CLASSIFY )
//#define BIT_WALK_METHOD scan_and_classify
//#elif  defined( USE_SCAN_AND_CLASSIFY_SWITCH )
//#define BIT_WALK_METHOD scan_and_classify_switch
//#elif  defined( USE_INLINE_AND_CLASSIFY )
//#define BIT_WALK_METHOD inline_classify
//#else   // default is to perform classification inline
//#define BIT_WALK_METHOD iterate_and_classify
//#endif  // classification procedure
//
//namespace clotho {
//namespace recombine {
//
//#define CHECK_0()   if( m_cfier( *first ) ) { res |= (Block) OFFSET( 0 ); }
//#define CHECK( x )  if( m_cfier( *(first + x) ) ) { res |= (Block) OFFSET( x ); }
//
//struct copy_matching_classify_mismatch {};
//
//template < class Classifier, class Tag = copy_matching_classify_mismatch >
//class bit_block_recombiner {
//public:
//    typedef Classifier classifier_type;
//
//    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}
//
//    template < class Block, class ElementIterator >
//    Block operator()( Block b0, Block b1, ElementIterator first ) {
//        return BIT_WALK_METHOD(b0, b1, first);
//    }
//
//    template < class Block, class ElementIterator >
//    Block iterate_and_classify( const Block b0, const Block b1, const ElementIterator first ) {
//        typedef clotho::utility::bit_block_iterator< Block, clotho::utility::tag::linear_iterator_tag > iterator;
//
//        Block mask = (Block)0;
//        iterator bit_it( (b0 ^ b1) ), bit_end;
//        while( bit_it != bit_end ) {
//            unsigned int idx = (*bit_it++);
//            if( m_cfier( *(first + idx) ) ) {
//                mask |= ((Block)1 << idx );
//            }
//        }
//
//        Block res = ((b0 & mask) | (b1 & ~ mask) );
//        return res;
//    }
//
//    /**
//     *  Inline classification walks each set bit in a block performing
//     *  the classification step inline
//     *
//     *  This results in:
//     *      - a single loop per 32-bit block (assuming 32 bit deBruijn hashing function)
//     *      - uses a switch to jump between hashed keys
//     *      - requires less bit shifting as mask offset bits can be inlined into case logic
//     *  Empirically performs efficiently, though slightly less efficiently than earlier linear bit block iterator
//     *  Still not as efficient as vector< index > approach.
//     *
//     *  Thought process:
//     *      - a 'single' iteration per 32-bit block is more efficient than multiple
//     *      - requires less memory as offsets are processed inline
//     *
//     *  Concerns:
//     *      - Inline analysis of bit positions may result in less streamline execution (more context switching)
//     */
//    template < class Block, class ElementIterator >
//    Block inline_classify( const Block b0, const Block b1, const ElementIterator first ) {
//        typedef Block block_type;
//        // Simplified Logic:
//        // b0 = 01001
//        // b1 = 10101
//        // (b0 ^ b1) = 11100
//        //
//        // assume:
//        // classify( bit_idx = 2 ) true
//        // classify( bit_idx = 3 ) true
//        // classify( bit_idx = 4 ) false
//        // => b0_mask = 01100
//        block_type b0_mask = (block_type) block_logic( (b0 ^ b1), first );
//
//        // b0 & b0_mask = (01001  &  01100) = 01000
//        // b1 & ~b0_mask = (10101 & ~01100) = (10101 & 10011) = 10001
//        // res = (01000 | 10001) = 11001
//        block_type res = ((b0 & b0_mask) | (b1 & ~b0_mask));
//
//        return res;
//    }
//
//    /**
//     * Scan and Classify walks the set bits in a block and stores their
//     * index (explodes the block into list of bit indices). The list of
//     * indices is iterated and the elements at these offsets are classified.
//     *
//     * This results in:
//     *      - a single bit walking loop per 32-bit block
//     *      - a simple loop over offsets
//     *      - switch over the bit indices
//     *
//     * Thought process:
//     *      - Bit scanning could more streamline (less context switching)
//     *      - Requires more memory (64 * 4 = 256 Bytes)
//     *      - Second loop could be more streamline
//     *      - Switch logic allows mask bits sets to be defined inline rather than computed
//     *
//     *  Concern:
//     *      - Doubles the number of iterative steps
//     */
//    template < class Block, class ElementIterator >
//    Block   scan_and_classify_switch( const Block b0, const Block b1, const ElementIterator first ) {
//        // Only compute the hash of each set bit
//        // In effect, delay the offset lookup to be performed by the switching logic
//        // Goal is to eliminate a 'double-lookup' (lookup from hash table, lookup based upon hash)
//        // This assumes switch logic results in a 'jump table'
//        unsigned int count = clotho::utility::hash_set_bits( (b0 ^ b1), indices );
//
//        Block res = (Block)0;
//        unsigned int * idx = indices;
//        while( count-- ) {
//            switch( *idx++ ) {
//            case 0:
//                CHECK_0() break;
//            case 1:
//                CHECK(1) break;
//            case 2:
//                CHECK(28) break;
//            case 3:
//                CHECK(2) break;
//            case 4:
//                CHECK(29) break;
//            case 5:
//                CHECK(14) break;
//            case 6:
//                CHECK(24) break;
//            case 7:
//                CHECK(3) break;
//            case 8:
//                CHECK(30) break;
//            case 9:
//                CHECK(22) break;
//            case 10:
//                CHECK(20) break;
//            case 11:
//                CHECK(15) break;
//            case 12:
//                CHECK(25) break;
//            case 13:
//                CHECK(17) break;
//            case 14:
//                CHECK(4) break;
//            case 15:
//                CHECK(8) break;
//            case 16:
//                CHECK(31) break;
//            case 17:
//                CHECK(27) break;
//            case 18:
//                CHECK(13) break;
//            case 19:
//                CHECK(23) break;
//            case 20:
//                CHECK(21) break;
//            case 21:
//                CHECK(19) break;
//            case 22:
//                CHECK(16) break;
//            case 23:
//                CHECK(7) break;
//            case 24:
//                CHECK(26) break;
//            case 25:
//                CHECK(12) break;
//            case 26:
//                CHECK(18) break;
//            case 27:
//                CHECK(6) break;
//            case 28:
//                CHECK(11) break;
//            case 29:
//                CHECK(5) break;
//            case 30:
//                CHECK(10) break;
//            case 31:
//                CHECK(9) break;
//            case 32:
//                CHECK(32) break;    // 0
//            case 33:
//                CHECK(33) break;    // 1
//            case 34:
//                CHECK(60) break;    // 28
//            case 35:
//                CHECK(34) break;    // 2
//            case 36:
//                CHECK(61) break;    // 29
//            case 37:
//                CHECK(46) break;    // 14
//            case 38:
//                CHECK(56) break;    // 24
//            case 39:
//                CHECK(35) break;    // 3
//            case 40:
//                CHECK(62) break;    // 30
//            case 41:
//                CHECK(54) break;    // 22
//            case 42:
//                CHECK(52) break;    // 20
//            case 43:
//                CHECK(47) break;    // 15
//            case 44:
//                CHECK(57) break;    // 25
//            case 45:
//                CHECK(49) break;    // 17
//            case 46:
//                CHECK(36) break;    // 4
//            case 47:
//                CHECK(40) break;    // 8
//            case 48:
//                CHECK(63) break;    // 31
//            case 49:
//                CHECK(59) break;    // 27
//            case 50:
//                CHECK(45) break;    // 13
//            case 51:
//                CHECK(55) break;    // 23
//            case 52:
//                CHECK(53) break;    // 21
//            case 53:
//                CHECK(51) break;    // 19
//            case 54:
//                CHECK(48) break;    // 16
//            case 55:
//                CHECK(39) break;    // 7
//            case 56:
//                CHECK(58) break;    // 26
//            case 57:
//                CHECK(44) break;    // 12
//            case 58:
//                CHECK(50) break;    // 18
//            case 59:
//                CHECK(38) break;    //  6
//            case 60:
//                CHECK(43) break;    // 11
//            case 61:
//                CHECK(37) break;    //  5
//            case 62:
//                CHECK(42) break;    // 10
//            case 63:
//                CHECK(41) break;    //  9
//            default:
//                break;
//            }
//        }
//
//        Block rec = ((b0 & res) | (b1 & ~res));
//        return rec;
//    }
//
///**
// *  Same logic as previous scan and classify
// *
// *  Thought process:
// *      - Simpler classification loop (indices will be ordered as a result of scan)
// *
// *  Concern:
// *      - Bit mask computation dependent upon index (adds bit shifting step)
// */
//    template < class Block, class ElementIterator >
//    Block scan_and_classify( const Block b0, const Block b1, const ElementIterator first ) {
//        unsigned int count = clotho::utility::scan_set_bits( (b0 ^ b1), indices );
//        unsigned int * idx = indices;
//
//        Block res = (Block)0;
//        while( count-- ) {
//            unsigned int tmp = *idx++;
//            if( m_cfier( *(first + tmp) ) ) {
//                res |= ((Block)1 << tmp);
//            }
//        }
//
//        Block rec = ((b0 & res) | (b1 & ~res));
//        return rec;
//    }
//
//    virtual ~bit_block_recombiner() {}
//
//protected:
//    classifier_type m_cfier;
//    unsigned int indices[ 64 ];
//
//    template < class ElementIterator >
//    inline unsigned long block_logic( unsigned long b, const ElementIterator first ) {
//        typedef unsigned long Block;
//
//        Block res = (Block)0;
//
//        if( !b ) return res;
//
//        unsigned int lo = (unsigned int) b;
//
//        while( lo ) {
//            unsigned int tmp = LEAST_SIG_BIT( lo );
//            switch( DEBRUIJNBIT_HASH( tmp ) ) {
//            case 0:
//                CHECK_0() break;
//            case 1:
//                CHECK(1) break;
//            case 2:
//                CHECK(28) break;
//            case 3:
//                CHECK(2) break;
//            case 4:
//                CHECK(29) break;
//            case 5:
//                CHECK(14) break;
//            case 6:
//                CHECK(24) break;
//            case 7:
//                CHECK(3) break;
//            case 8:
//                CHECK(30) break;
//            case 9:
//                CHECK(22) break;
//            case 10:
//                CHECK(20) break;
//            case 11:
//                CHECK(15) break;
//            case 12:
//                CHECK(25) break;
//            case 13:
//                CHECK(17) break;
//            case 14:
//                CHECK(4) break;
//            case 15:
//                CHECK(8) break;
//            case 16:
//                CHECK(31) break;
//            case 17:
//                CHECK(27) break;
//            case 18:
//                CHECK(13) break;
//            case 19:
//                CHECK(23) break;
//            case 20:
//                CHECK(21) break;
//            case 21:
//                CHECK(19) break;
//            case 22:
//                CHECK(16) break;
//            case 23:
//                CHECK(7) break;
//            case 24:
//                CHECK(26) break;
//            case 25:
//                CHECK(12) break;
//            case 26:
//                CHECK(18) break;
//            case 27:
//                CHECK(6) break;
//            case 28:
//                CHECK(11) break;
//            case 29:
//                CHECK(5) break;
//            case 30:
//                CHECK(10) break;
//            case 31:
//                CHECK(9) break;
//            default:
//                break;
//            }
//            lo ^= tmp;
//        }
//
//        lo = (unsigned int) (b >> 32);
//        while( lo ) {
//            unsigned int tmp = LEAST_SIG_BIT( lo );
//            switch( DEBRUIJNBIT_HASH( tmp ) ) {
//            case 0:
//                CHECK(32) break;    // 0
//            case 1:
//                CHECK(33) break;    // 1
//            case 2:
//                CHECK(60) break;    // 28
//            case 3:
//                CHECK(34) break;    // 2
//            case 4:
//                CHECK(61) break;    // 29
//            case 5:
//                CHECK(46) break;    // 14
//            case 6:
//                CHECK(56) break;    // 24
//            case 7:
//                CHECK(35) break;    // 3
//            case 8:
//                CHECK(62) break;    // 30
//            case 9:
//                CHECK(54) break;    // 22
//            case 10:
//                CHECK(52) break;    // 20
//            case 11:
//                CHECK(47) break;    // 15
//            case 12:
//                CHECK(57) break;    // 25
//            case 13:
//                CHECK(49) break;    // 17
//            case 14:
//                CHECK(36) break;    // 4
//            case 15:
//                CHECK(40) break;    // 8
//            case 16:
//                CHECK(63) break;    // 31
//            case 17:
//                CHECK(59) break;    // 27
//            case 18:
//                CHECK(45) break;    // 13
//            case 19:
//                CHECK(55) break;    // 23
//            case 20:
//                CHECK(53) break;    // 21
//            case 21:
//                CHECK(51) break;    // 19
//            case 22:
//                CHECK(48) break;    // 16
//            case 23:
//                CHECK(39) break;    // 7
//            case 24:
//                CHECK(58) break;    // 26
//            case 25:
//                CHECK(44) break;    // 12
//            case 26:
//                CHECK(50) break;    // 18
//            case 27:
//                CHECK(38) break;    //  6
//            case 28:
//                CHECK(43) break;    // 11
//            case 29:
//                CHECK(37) break;    //  5
//            case 30:
//                CHECK(42) break;    // 10
//            case 31:
//                CHECK(41) break;    //  9
//            default:
//                break;
//            }
//            lo ^= tmp;
//        }
//        return res;
//    }
//};
//
//#undef CHECK_0
//#undef CHECK
//
//}   // namespace recombinations {
//}   // namespace clotho {

#endif  // BLOCK_RECOMBINER_HPP_
