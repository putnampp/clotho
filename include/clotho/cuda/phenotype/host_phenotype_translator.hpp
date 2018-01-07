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
        if( traits.isAllNeutral() ) {
            evaluate_phenotype::reset( pop->getDevicePhenotypes(), pop->getPhenotypeCapacity() );
        } else {
            typename evaluate_phenotype::warp_sequence_kernel v;
            evaluate_phenotype::execute( pop->getDeviceSequences(), traits.getDeviceWeights(), pop->getDevicePhenotypes(), pop->getSequenceCount(), pop->getBlocksPerSequence(), all_count, traits.getAlleleCount(), traits.getTraitCount(), &v );
        }
    }

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop, HostTraitSpace< RealType > & traits, unsigned int all_count, cudaStream_t & stream ) {
        if( traits.isAllNeutral() ) {
            evaluate_phenotype::reset( pop->getDevicePhenotypes(), pop->getPhenotypeCapacity() );
        } else {
            typename evaluate_phenotype::warp_sequence_kernel v;
            evaluate_phenotype::execute( pop->getDeviceSequences(), traits.getDeviceWeights(), pop->getDevicePhenotypes(), pop->getSequenceCount(), pop->getBlocksPerSequence(), all_count, traits.getAlleleCount(), traits.getTraitCount(), stream, &v );
        }
    }

    virtual ~HostPhenotypeTranslator() {}
};

#endif  // HOST_PHENOTYPE_TRANSLATOR_HPP_
