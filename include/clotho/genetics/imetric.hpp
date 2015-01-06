#ifndef IMETRIC_HPP_
#define IMETRIC_HPP_

#include <string>

class ifitness {
public:

    virtual double operator()( double ) = 0;
    virtual const std::string name() const = 0;
};

#endif  // IMETRIC_HPP_
