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
#ifndef ALLELE_FREQUENCY_HPP_
#define ALLELE_FREQUENCY_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/cuda/data_spaces/basic_data_space.hpp"
#include "clotho/cuda/analysis/allele_frequency_kernels.hpp"

#include "clotho/cuda/device_state_object.hpp"

class AlleleFrequency : public clotho::utility::iStateObject {
public:

    AlleleFrequency( boost::property_tree::ptree & config ) :
        blocks(200, 1, 1)
        , threads( 32, 32, 1)
    {
        create_space( dCounts );
        parse_configuration( config );
    }

    template < class PopulationSpaceType >
    void evaluate( PopulationSpaceType * pop ) {
        _resize_space<<< 1, 1 >>>( dCounts, pop->alleles.get_device_space() );
        count_alleles<<< blocks, threads >>>( pop->sequences.get_device_space(), dCounts ); 
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree f;

        get_device_object_state( f, dCounts );

        state.add_child( "allele_frequency", f );
    }

    virtual ~AlleleFrequency() {
        delete_space( dCounts );
    }

protected:
    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lc, bconf, tconf;

        lc = config.get_child( "allele_frequency", lc );

        bconf.put( "x", blocks.x );
        bconf.put( "y", blocks.y );

        tconf.put( "x", threads.x );
        tconf.put( "y", threads.y );
        
        lc.put_child( "kernel.blocks", bconf );
        lc.put_child( "kernel.threads", tconf );

        config.put_child( "allele_frequency", lc );
    }

    basic_data_space< unsigned int > * dCounts;
    dim3 blocks, threads;
};

#endif  // ALLELE_FREQUENCY_HPP_
