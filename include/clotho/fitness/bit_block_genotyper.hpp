#ifndef CLOTHO_BIT_BLOCK_GENOTYPER_HPP_
#define CLOTHO_BIT_BLOCK_GENOTYPER_HPP_

#include <type_traits>

namespace clotho {
namespace fitness {

template < class Block, class Enabled = void >
class bit_block_genotyper;

template < class Block >
class bit_block_genotyper< Block, typename std::enable_if< std::is_integral< Block >::value >::type >  {
public:
    typedef Block block_type;

    static block_type   get_all_heterozygous( const block_type b0, const block_type b1 ) {
        return (b0 ^ b1);
    }

    static block_type get_all_homozygous( const block_type b0, const block_type b1 ) {
        return (~(b0 ^ b1));
    }

    static block_type   get_alt_homozygous( const block_type b0, const block_type b1 ) {
        return (b0 & b1);
    }

    static block_type   get_ref_homozygous( const block_type b0, const block_type b1 ) {
        return (~(b0 | b1));
    }
};

}   // namespace fitness {
}   // namespace clotho {

#endif  // CLOTHO_BIT_BLOCK_GENOTYPER_HPP_
