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

#ifndef CUDA_PHENOTYPE_TRANSLATE_ALGORITHM
#define CUDA_PHENOTYPE_TRANSLATE_ALGORITHM 2
#endif

class phenotype_translator {
public:

    phenotype_translator( boost::property_tree::ptree & config ) {
        parse_configuration( config );
    }

    template < class PopulationSpaceType >
    void translate( PopulationSpaceType * pop, unsigned int seq_count ) {
        // algorithm 2 performs better than algorithm 3
        // algorithm 3 performs better than algorithm 1
        //
        clotho::utility::algo_version< CUDA_PHENOTYPE_TRANSLATE_ALGORITHM > * v = NULL;
        translate_algo( pop, seq_count, v );
    }

    virtual ~phenotype_translator() {}
protected:

    void parse_configuration( boost::property_tree::ptree & config ) {

    }

    template < class PopulationSpaceType, unsigned char V >
    void translate_algo( PopulationSpaceType * pop, unsigned int seq_count, clotho::utility::algo_version< V > * v) {
        // blocks => x == sequences; y == phenotype traits
        // threads => x == alleles per sequence block; y == parallel sequence blocks to operate on 
        dim3 blocks( 100, 1, 1 ), threads( 32, 32, 1 );
        _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v );
    }

    template < class PopulationSpaceType >
    void translate_algo( PopulationSpaceType * pop, unsigned int seq_count, clotho::utility::algo_version< 1 > * v ) {
        dim3 blocks( seq_count, 1, 1), threads( 32, 32, 1 );
        _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v, 0, seq_count );
    }

    template < class PopulationSpaceType >
    void translate_algo( PopulationSpaceType * pop, unsigned int seq_count, clotho::utility::algo_version< 2 > * v ) {
        unsigned int b = seq_count / 32;    // full blocks (sequences)

        dim3 blocks( b, 1, 1 ), threads( 32, 32, 1 );
        _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v, 0, b * 32 );

        if( seq_count % 32 != 0 ) {
            dim3 tail_block( 1, 1, 1), tail_warps( 32, seq_count % 32, 1 );

            _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v, b * 32, seq_count );
        }
    }

    template < class PopulationSpaceType >
    void translate_algo( PopulationSpaceType * pop, unsigned int seq_count, clotho::utility::algo_version< 3 > * v ) {
        unsigned int b = (seq_count / (32 * 32)) + ((seq_count % (32 * 32) == 0) ? 0 : 1);

        dim3 blocks( b, 1, 1), threads( 32, 32, 1 );
        _translate<<< blocks, threads >>>( pop->alleles.get_device_space(), pop->sequences.get_device_space(), pop->pheno_space, v, 0, seq_count );
    }
};
#endif  // PHENOTYPE_TRANSLATOR_HPP_
