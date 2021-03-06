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
#ifndef DEVICE_PHENOTYPE_SPACE_DEF_HPP_
#define DEVICE_PHENOTYPE_SPACE_DEF_HPP_

/*
template < class RealType >
struct device_phenotype_space {
    typedef RealType    real_type;

    real_type   * phenotypes;

    unsigned int size, capacity;
};*/


#include "clotho/cuda/data_spaces/basic_data_space.hpp"

template < class RealType >
struct device_phenotype_space : public basic_data_space< RealType > {
    typedef typename basic_data_space< RealType >::value_type    real_type;
};

#endif  // DEVICE_PHENOTYPE_SPACE_DEF_HPP_
