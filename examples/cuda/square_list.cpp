#include <stdio.h>
#include <cuda.h>
#include <iostream>

#include <boost/random/mersenne_twister.hpp>

#include "clotho/utility/timer.hpp"

#include "square.h"

typedef clotho::utility::timer  timer_type;
typedef boost::random::mt19937  rng_type;

typedef Square::int_type        int_type;

void square_it( int_type N ) {

    timer_type t;
    int_type * a = new int_type[ N ];

    for(int_type i = 0; i < N; ++i ) {
        a[i] = i*i;
    }

    delete [] a;
    t.stop();
    std::cout << "host lapse:   " << t.elapsed().count() << std::endl;
}

void square_list_rand( int_type N ) {
    timer_type t;

    int_type * a = new int_type[ N ];
    rng_type rand(1234);

    for(int_type i = 0; i < N; ++i ) {
        unsigned int r = rand();

        a[i] = r * r;
    }

    t.stop();
    std::cout << "host random lapse: " << t.elapsed().count() << std::endl;

    std::cerr << "BEGIN HOST" << std::endl;
    for( unsigned int i = 0; i < N; ++i ) {
        std::cerr << i << " -> " << a[i] << std::endl;
    }
    delete [] a;
}

int main( int argc, char ** argv ) {

    timer_type t;
    Square s;
    t.stop();

    std::cout << "Init: " << t.elapsed().count() << std::endl;

    //t.start();
    //s();
    //t.stop();

    //std::cout << "device lapse: " << t.elapsed().count() << std::endl;

    t.start();
    s.random_list();
    t.stop();
    std::cout << "device random lapse: " << t.elapsed().count() << std::endl;
    std::cerr << s;

    square_it( s.size() );

    square_list_rand( s.size() );

    return 0;
}
