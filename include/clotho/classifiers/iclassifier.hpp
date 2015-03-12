#ifndef ICLASSIFIER_HPP_
#define ICLASSIFIER_HPP_

#include <string>

template < class Element >
struct iclassifier {

    virtual const std::string name() const = 0;
    virtual bool operator()( const Element & elem ) const = 0;

    virtual ~iclassifier();
};

#endif  // ICLASSIFIER_HPP_
