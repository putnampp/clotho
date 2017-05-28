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
        dim3 blocks( 1,1,1), threads( 1,1,1 );
        computBlockThreadDims( pop->getBlocksPerSequence(), blocks, threads );

        evaluate_free_space<<< blocks, threads >>>( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence() );
    }

    template < class RealType, class IntType >
    void remove_fixed( HostPopulationSpace< RealType, IntType > * pop ) {
        dim3 blocks( 1,1,1), threads( 1,1,1 );
        computBlockThreadDims( pop->getBlocksPerSequence(), blocks, threads );

        remove_fixed<<< blocks, threads >>>( pop->getDeviceSequences(), pop->getDeviceFreeSpace(), pop->getSequenceCount(), pop->getBlocksPerSequence() );
    }

    virtual ~HostFreeSpace() {}

protected:
    void computeBlockThreadDims( unsigned int N, dim3 & blocks, dim3 & threads ) {
        assert( N % 32 == 0);

        unsigned int t_y = N / 32;
        unsigned int b_x = 1;

        if( t_y > 32 ) {
            b_x = t_y / 32;
            if( t_y % 32 != 0 ) {
                b_x += 1;
            }
            t_y = 32;
        }

        blocks.x = b_x;
        threads.y = t_y;
    }
};
#endif  // HOST_FREE_SPACE_HPP_
