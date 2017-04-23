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
#ifndef PHENOTYPE_TRANSLATOR_HPP_
#define PHENOTYPE_TRANSLATOR_HPP_

#include <cuda.h>
#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/phenotype/translation_kernels.hpp"

class phenotype_translator {
public:

    phenotype_translator( boost::property_tree::ptree & config ) {
        parse_configuration( config );
    }

    template < class PopulationSpaceType >
    void translate( PopulationSpaceType * pop ) {
        clotho::utility::algo_version< 1 > * v = NULL;

        // blocks => x == sequences; y == phenotype traits
        // threads => x == alleles per sequence block; y == parallel sequence blocks to operate on 
        dim3 blocks( 100, 1, 1 ), threads( 32, 32, 1 );
        _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v );
    }

    virtual ~phenotype_translator() {}
protected:

    void parse_configuration( boost::property_tree::ptree & config ) {

    }
};
#endif  // PHENOTYPE_TRANSLATOR_HPP_
