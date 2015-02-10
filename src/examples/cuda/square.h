#ifndef SQUARE_CUDA_H_
#define SQUARE_CUDA_H_

#include <ostream>

class Square {
public:
    static const unsigned int N = 15;

    Square();
    void operator()();

    friend std::ostream & operator<<( std::ostream &, const Square & rhs );

    virtual ~Square();
protected:
    unsigned int m_a[N];
};

#endif  // SQUARE_CUDA_H_
