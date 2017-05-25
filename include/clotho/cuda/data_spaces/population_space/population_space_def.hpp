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
#include "clotho/cuda/data_spaces/phenotype_space/device_phenotype_space.hpp"

#include "clotho/cuda/data_spaces/free_space/device_free_space.hpp"
#include "clotho/cuda/data_spaces/event_space/device_event_space.hpp"

#include "clotho/cuda/data_spaces/population_space/population_space_helper_api.hpp"
#include "clotho/cuda/device_state_object.hpp"

template < class RealType, class IntType, class OrderTag >
struct PopulationSpace : public clotho::utility::iStateObject {

    typedef PopulationSpace< RealType, IntType, OrderTag >  self_type;

    typedef SequenceSpace< IntType >                        sequence_space_type;
    typedef AlleleSpace< RealType >                         allele_space_type;

//    typedef device_phenotype_space< RealType >              phenotype_space_type;
    typedef basic_data_space< RealType >                    phenotype_space_type;
    typedef device_free_space< IntType, OrderTag >          free_space_type;
    typedef device_event_space< IntType, OrderTag >         event_space_type;

    typedef typename allele_space_type::real_type   real_type;
    typedef typename sequence_space_type::int_type  int_type;
    typedef typename free_space_type::order_tag_type         order_tag_type;
    
    sequence_space_type sequences;
    allele_space_type   alleles;

    free_space_type         * free_space;
    phenotype_space_type    * pheno_space;

    unsigned int            seq_count;

    PopulationSpace( boost::property_tree::ptree & config ) :
        sequences( config )
        , alleles( config )
        , seq_count(0)
        {
        create_space( free_space );
        create_space( pheno_space );

//        resize_space( free_space, allele_space_type::device_space_type::ALIGNMENT_SIZE );

        update_metadata();
    }

    void resize( self_type * parent_pop, event_space_type * mut_events, unsigned int seqs ) {
        resize_space( pheno_space, seqs); 

        alleles.expand_relative_to( parent_pop->alleles, parent_pop->free_space, free_space, mut_events );
        cudaDeviceSynchronize();

        update_free_map( free_space );

        validate_free_map( free_space, mut_events );

        sequences.resize( alleles.get_device_space(), seqs );

        seq_count = seqs;
    }

    void update_metadata() {
        update_free_space2( sequences.get_device_space(), free_space );
    }

    void record_fixed( allele_space_type & fixed ) {
        alleles.move_fixed( fixed, free_space );
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree a, s, f, p;

        alleles.get_state( a );
        sequences.get_state( s );

        get_device_object_state( f, free_space );
        get_device_object_state( p, pheno_space );

        state.put_child( "alleles", a );
        state.put_child( "sequences", s );
        state.put_child( "free_space", f );
        state.put_child( "phenotypes", p );
    }

    virtual ~PopulationSpace() {
        delete_space( free_space );
        delete_space( pheno_space );
    }
};

#endif  // POPULATION_SPACE_DEF_HPP_
