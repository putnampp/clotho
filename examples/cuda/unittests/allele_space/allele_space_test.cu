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
#include <iostream>
#include <cassert>

#include "clotho/cuda/allele_space/allele_space.hpp"
#include "clotho/utility/timer.hpp"

typedef double          real_type;
typedef unsigned int    int_type;

#ifdef USE_UNIT_ORDERED
typedef unit_ordered_tag< int_type > tag_type;
#else
typedef unordered_tag tag_type;
#endif  // USE_UNIT_ORDERED

typedef clotho::utility::timer timer_type;

int main( int argc, char ** argv ) {
    const unsigned int COUNT = 3;
    unsigned int sizes[COUNT] = { 128, 64, 1024 };

    for( unsigned int i = 0; i < COUNT; ++i ) {
        AlleleSpace< real_type, int_type, tag_type > as;
        timer_type t;
        as.resize( sizes[i] );
        t.stop();
        std::cout << "Resize lag: " << t << std::endl;

        std::cout << as << std::endl;
    }

    return 0;
}
