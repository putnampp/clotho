//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
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
