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

        _translate<<< 200, 128 >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space );
    }

    virtual ~phenotype_translator() {}
protected:

    void parse_configuration( boost::property_tree::ptree & config ) {

    }
};
#endif  // PHENOTYPE_TRANSLATOR_HPP_
