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
#ifndef HOST_FREE_SPACE_HPP_
#define HOST_FREE_SPACE_HPP_

#include "clotho/cuda/free_space/free_space_kernel.hpp"
#include "clotho/cuda/analysis/sequence_weight_kernel.hpp"

class HostFreeSpace {
public:

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop, unsigned int max_index ) {
        evaluate_free_space::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence());

        dim3 blocks(3, 1, 1), threads( 32, 32, 1 );
        evaluate_sequence_weights<<< blocks, threads >>>( pop->getDeviceFreeSpace(), pop->getDeviceCounts(), pop->getBlocksPerSequence(), max_index ); 
    }

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop, unsigned int max_index, cudaStream_t & stream ) {
        evaluate_free_space::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence(), stream);

        dim3 blocks(3, 1, 1), threads( 32, 32, 1 );
        evaluate_sequence_weights<<< blocks, threads, 0, stream >>>( pop->getDeviceFreeSpace(), pop->getDeviceCounts(), pop->getBlocksPerSequence(), max_index ); 
    }

    template < class RealType, class IntType >
    void perform_remove_fixed( HostPopulationSpace< RealType, IntType > * pop ) {
        remove_fixed::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence() );
    }

    virtual ~HostFreeSpace() {}

protected:
};
#endif  // HOST_FREE_SPACE_HPP_
