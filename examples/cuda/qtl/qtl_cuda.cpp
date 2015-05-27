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

#include "device_properties.h"
#include "page_manager.h"

#include "adjacency_matrix.h"

int main(int argc, char ** argv ) {

    std::cout << "Hello world" << std::endl;

    DeviceProperties props;

    std::cout << props << std::endl;   

    adjacency_matrix mat( 225, 32 );

    std::cout << mat << std::endl;

    assert( mat.capacity() == 2 );

    mat.resize( 4000, 5000 );

    std::cout << mat;

    return 0;
}
