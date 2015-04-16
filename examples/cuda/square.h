#ifndef SQUARE_CUDA_H_
#define SQUARE_CUDA_H_

#include <ostream>

class Square {
public:
    typedef unsigned long int_type;

    Square();

    size_t size() const;

    void operator()();

    void random_list();

    friend std::ostream & operator<<( std::ostream &, const Square & rhs );

    virtual ~Square();
protected:

    void init();
    //int_type m_a[N];

    int_type *  m_a, * m_dest;
    size_t      m_size;
    int m_maxBlocks, m_maxThreadsPerBlock;
};

#endif  // SQUARE_CUDA_H_
