#include <stdio.h>
#include <cuda.h>
#include <iostream>

//#include "clotho/utility/timer.hpp"

#include "square.h"

//typedef clotho::utility::timer  timer_type;

int main( int argc, char ** argv ) {

//    unsigned int * a_h, * a_d;
//
//    size_t size = N * sizeof(unsigned int);
//    a = (unsigned int *) malloc(size);
//    a_h = new unsigned int[N];
//
//    timer_type t;
//
//    cudaMalloc( (void **) &a_d, size);
//
//    t.stop();
//
//    std::cout << "Malloc: " << t.elasped().count() << std::endl;
// 
//    t.start();
//    square<<<1,N>>>(a_d, N);
//    t.stop();
//
//    std::cout << "Square: " << t.elapsed().count() << std::endl;
//
//    cudaMemcpy(a_h, a_d, size, cudaMemcpyDeviceToHost);

//    for( unsigned int i = 0; i < N; ++i ) {
//        std::cout << i << " -> " << a_h[i] << std::endl;
//    }
//
//    cudaFree(a_d);
//  free( a_h );
//    delete [] a_h;

    Square s;
    s();
    std::cout << s;

    return 0;
}
