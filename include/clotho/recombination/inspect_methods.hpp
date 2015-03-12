#ifndef INSPECT_METHODS_HPP_
#define INSPECT_METHODS_HPP_

namespace clotho {
namespace recombine {
namespace inspection {
namespace tag {

struct copy_matching_classify_mismatch {
    template < class Block >
    static Block select( const Block b0, const Block b1 ) {
        return (b0 ^ b1);
    }
};

struct classify_all {
    template < class Block >
    static Block select( const Block b0, const Block b1 ) {
        return (b0 | b1);
    }
};

}   // namespace tag
}   // namespace inspection
}   // namespace recombine
}   // namespace clotho

#endif  // INSPECT_METHODS_HPP_
