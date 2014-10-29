#ifndef BIT_BLOCK_RECOMBINER_HPP_
#define BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/utility/bit_block_iterator.hpp"

namespace clotho {
namespace recombine {

struct copy_matching_classify_mismatch {};

template < class Classifier, class Tag = copy_matching_classify_mismatch >
class bit_block_recombiner {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( Block b0, Block b1, ElementIterator first ) {
        typedef Block block_type;
        typedef clotho::utility::bit_block_iterator< Block > iterator;

        Block matching = (b0 & b1);
        Block mismatch = (b0 ^ b1);

        Block b0_mask = (Block)0;
        iterator bit_it( mismatch ), bit_end;
        while( bit_it != bit_end ) {
            unsigned int bit_idx = (*bit_it++);
            if( m_cfier(first, bit_idx) ) {
                b0_mask |= ((block_type)1 << (bit_idx));
            }
        }

        Block b1_mismatch = ((~b0_mask & mismatch) & b1);
        Block res = matching | b1_mismatch;
        res |= (b0_mask & b0);
        return  res;
    }

    virtual ~bit_block_recombiner() {}

protected:
    classifier_type m_cfier;

};

}   // namespace recombinations {
}   // namespace clotho {

#endif  // BLOCK_RECOMBINER_HPP_
