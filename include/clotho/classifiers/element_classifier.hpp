#ifndef ELEMENT_CLASSIFIER_HPP_
#define ELEMENT_CLASSIFIER_HPP_

namespace clotho {
namespace classifier {


template < class Result = bool >
class element_classifier {
public:
    typedef Result  result_type;

    template < class Element >
    result_type operator()( const element_type & elem ) const {
        return false;
    }

    template < class ElementIterator >
    result_type operator()( ElementIterator elem_it, size_t idx) const {
        return operator()( *(elem_it + idx) );
    }
};

}   // namespace classifier {
}   // namespace clotho {

#endif  // ELEMENT_CLASSIFIER_HPP_
