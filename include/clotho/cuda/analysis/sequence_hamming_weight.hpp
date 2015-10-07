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
#ifndef SEQUENCE_HAMMING_WEIGHT_HPP_
#define SEQUENCE_HAMMING_WEIGHT_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include "clotho/cuda/analysis/sequence_hamming_weight_kernel.hpp"

#include "clotho/cuda/device_state_object.hpp"

#ifndef HAMMING_VERSION
#define HAMMING_VERSION 4
#endif  // HAMMING_VERSION

class SequenceHammingWeight : public clotho::utility::iStateObject {
public:

    typedef basic_data_space< unsigned int >    space_type;
    typedef algo_version< HAMMING_VERSION >    algo_version_type;

    SequenceHammingWeight( boost::property_tree::ptree & config ) :
        dWeights( NULL )
        , blocks( 200, 1, 1 )
        , threads( 32, 32, 1)
    {
        create_space( dWeights );
        parse_configuration( config );
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop ) {
        _resize_space_to<<< 1, 1 >>>( dWeights, pop->sequences.get_device_space() );
        CHECK_LAST_KERNEL_EXEC

        algo_version_type * v = NULL;

        sequence_hamming_weight_kernel<<< blocks, threads >>>( pop->sequences.get_device_space(), dWeights, v );
        CHECK_LAST_KERNEL_EXEC
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree w;
        get_device_object_state( w, dWeights );
        state.add_child( "sequence_degree_spectrum", w );
    }

    virtual ~SequenceHammingWeight() {
        delete_space( dWeights );
    }

protected:
    void parse_configuration( boost::property_tree::ptree & config ) {

    }

    space_type * dWeights;
    dim3 blocks, threads;
};

#endif  // SEQUENCE_HAMMING_WEIGHT_HPP_
