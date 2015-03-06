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
struct inline_dynamic_classify {};

}   // namespace tag
}   // namespace walker
}   // namespace recombine
}   // namespace clotho

namespace clotho {
namespace recombine {

#define CHECK_0()   if( m_cfier( *first ) ) { res |= (Block) OFFSET( 0 ); }
#define CHECK( x )  if( m_cfier( *(first + x) ) ) { res |= (Block) OFFSET( x ); }

template < class Classifier, class InspectMethodTag >
class bit_block_recombiner< Classifier, InspectMethodTag, clotho::recombine::walker::tag::inline_dynamic_classify > {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( const Block b0, const Block b1, const ElementIterator first ) {
        Block b = InspectMethodTag::select(b0, b1);

        Block mask = walk( b, first );

        Block rec = ((b0 & mask) | (b1 & ~mask));

        return rec;
    }

protected:
    template < class ElementIterator >
    unsigned int walk( unsigned int b, const ElementIterator first ) {
        unsigned int res = 0;
        while( b ) {
            unsigned int tmp = LEAST_SIG_BIT( b );
            unsigned int offset = DEBRUIJNBIT_HASH_LOOKUP( tmp );
            if( m_cfier( *(first + offset) ) ) {
                res |= ((Block)1 << offset);
            }
            b ^= tmp;
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

    classifier_type m_cfier;
};

#undef CHECK_0
#undef CHECK

}   // namespace recombine
}   // namespace clotho
#endif  // INLINE_BIT_BLOCK_RECOMBINER_HPP_
