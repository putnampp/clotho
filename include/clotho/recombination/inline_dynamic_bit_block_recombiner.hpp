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

    template < class ElementIterator >
    unsigned int operator( const unsigned int b0, const unsigned int b1, const ElementIterator first ) {
        typedef unsigned int Block;

        Block res = (Block)0;

        Block b = InspectMethodTag::select(b0, b1);

        while( b ) {
            unsigned int tmp = LEAST_SIG_BIT( b );
            unsigned int offset = DEBRUIJNBIT_HASH_LOOKUP( tmp );
            if( m_cfier( *(first + offset) ) ) {
                res |= ((Block)1 << offset);
            }
            b ^= tmp;
        }
        Block rec = ((b0 & res) | (b1 & ~res));

        return rec;
    }

    template < class ElementIterator >
    unsigned long operator()( const unsigned long b0, const unsigned long b1, const ElementIterator first ) {
        typedef unsigned long Block;

        Block res = (Block)0;

        Block b = InspectMethodTag::select(b0, b1);
        if( b ) {
            unsigned int lo = (unsigned int) b;

            while( lo ) {
                unsigned int tmp = LEAST_SIG_BIT( lo );
                unsigned int offset = DEBRUIJNBIT_HASH_LOOKUP( tmp );
                if( m_cfier( *(first + offset) ) ) {
                    res |= ((Block)1 << offset);
                }
                lo ^= tmp;
            }

            lo = (unsigned int) (b >> 32);
            Block hi_mask = (Block)0;
            while( lo ) {
                unsigned int tmp = LEAST_SIG_BIT( lo );
                unsigned int offset = DEBRUIJNBIT_HASH_LOOKUP( tmp );
                if( m_cfier( *(first + offset) ) ) {
                    hi_mask |= ((Block)1 << offset);
                }
                lo ^= tmp;
            }

            res |= (hi_mask << 32);
        }

        Block rec = ((b0 & res) | (b1 & ~res));

        return rec;
    }

protected:
    classifier_type m_cfier;
};

#undef CHECK_0
#undef CHECK

}   // namespace recombine
}   // namespace clotho
#endif  // INLINE_BIT_BLOCK_RECOMBINER_HPP_
