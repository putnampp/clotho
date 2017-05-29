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
#ifndef HOST_PHENOTYPE_TRANSLATOR_HPP_
#define HOST_PHENOTYPE_TRANSLATOR_HPP_

#include "clotho/cuda/phenotype/phenotype_evaluator.hpp"

class HostPhenotypeTranslator {
public:

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop, HostTraitSpace< RealType > & traits, unsigned int all_count ) {
        traits.updateDevice();

        assert( traits.getAlleleCount() % 32 == 0 );
        assert( pop->getDeviceSequences() != NULL );
        assert( traits.getDeviceWeights() != NULL );
        assert( pop->getDevicePhenotypes() != NULL );

        dim3 blocks( pop->getSequenceCount(), traits.getTraitCount(), 1 ), threads( 32, 32, 1 );
        evaluate_phenotype<<< blocks, threads >>>( pop->getDeviceSequences(), traits.getDeviceWeights(), pop->getDevicePhenotypes(), pop->getBlocksPerSequence(), all_count, traits.getAlleleCount() );
    }

    virtual ~HostPhenotypeTranslator() {}
};

#endif  // HOST_PHENOTYPE_TRANSLATOR_HPP_
