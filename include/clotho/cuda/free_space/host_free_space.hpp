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

class HostFreeSpace {
public:

    template < class RealType, class IntType >
    void operator()( HostPopulationSpace< RealType, IntType > * pop ) {
        evaluate_free_space::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence());
    }

    void operator()( HostPopulationSpace< RealType, IntType > * pop, cudaStream_t & stream ) {
        evaluate_free_space::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence(), stream);
    }

    template < class RealType, class IntType >
    void perform_remove_fixed( HostPopulationSpace< RealType, IntType > * pop ) {
        remove_fixed::execute( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence() );
    }

    virtual ~HostFreeSpace() {}

protected:
};
#endif  // HOST_FREE_SPACE_HPP_
