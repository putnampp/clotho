#ifndef FITNESS_METHOD_HPP_
#define FITNESS_METHOD_HPP_

namespace clotho {
namespace fitness {

struct multiplicative_heterozygous_tag {};
struct multiplicative_homozygous_tag {};

struct additive_heterozygous_tag {};
struct additive_homozygous_tag {};

template < class Result, class Tag >
class fitness_method {
public:
    typedef Result result_type;
    
    template < class Element >
    void operator()( result_type & res, const Element & elem, result_type scale = 1. ) {}
};

}   // namespace fitness {
}   // namespace clotho {

#endif  // FITNESS_METHOD_HPP_
