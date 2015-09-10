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
#ifndef POPULATION_SPACE_DEF_HPP_
#define POPULATION_SPACE_DEF_HPP_

#include "clotho/cuda/data_spaces/allele_space/allele_space.hpp"
#include "clotho/cuda/data_spaces/sequence_space/sequence_space.hpp"
#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"

#include "clotho/cuda/data_spaces/population_space/population_space_helper_api.hpp"

template < class RealType, class IntType, class OrderTag >
struct PopulationSpace {

    typedef PopulationSpace< RealType, IntType, OrderTag >  self_type;

    typedef SequenceSpace< IntType >                        sequence_space_type;
    typedef AlleleSpace< RealType, IntType, OrderTag >      allele_space_type;
    typedef device_free_space< IntType, OrderTag >          free_space_type;
    
    sequence_space_type sequences;
    allele_space_type   alleles;

    free_space_type    * free_space;

    PopulationSpace() {
        create_space( free_space );
        sequences.resize( alleles.get_device_space(), 2 );
    }

    template < class EventSpaceType >
    void resize( self_type * parent_pop, EventSpaceType * mut_events, unsigned int seqs ) {

//        std::cerr << "Resizing population" << std::endl;
        update_free_space( parent_pop->sequences.get_device_space(), free_space );

        alleles.expand_relative_to( parent_pop->alleles, free_space, mut_events );

        sequences.resize( alleles.get_device_space(), seqs );
    }

    virtual ~PopulationSpace() {
        delete_space( free_space );
    }
};

#endif  // POPULATION_SPACE_DEF_HPP_
