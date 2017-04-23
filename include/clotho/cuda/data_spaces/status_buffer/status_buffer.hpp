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
#ifndef CLOTHO_STATUS_BUFFER_HPP_
#define CLOTHO_STATUS_BUFFER_HPP_

#include "clotho/cuda/device_state_object.hpp"
#include "clotho/cuda/data_spaces/basic_data_space.hpp"

#include <boost/property_tree/ptree.hpp>

template < class PopulationSpaceType >
class StatusBuffer : public clotho::utility::iStateObject {
public:
    typedef PopulationSpaceType                 population_type;
    typedef basic_data_space< unsigned int >    count_buffer_type;

    StatusBuffer( unsigned int N ) {
        create_space( free_space_counts );
        create_space( allele_space_sizes );

        resize_space( free_space_counts, N );
        resize_space( allele_space_sizes, N );
    }

    void track( population_type * pop, unsigned int idx ) {
        record_size<<< 1, 1 >>>( pop->free_space, free_space_counts, idx );
        record_size<<< 1, 1 >>>( pop->alleles.get_device_space(), allele_space_sizes, idx );
    }

    void get_state( boost::property_tree::ptree & state ) {
        boost::property_tree::ptree fc, ac;

        get_device_object_state( fc, free_space_counts );
        get_device_object_state( ac, allele_space_sizes );

        state.put_child( "allele_count", ac );
        state.put_child( "free_count", fc );

    }

    ~StatusBuffer() {
        delete_space( free_space_counts );
        delete_space( allele_space_sizes );
    }

protected:
    count_buffer_type   * free_space_counts;
    count_buffer_type   * allele_space_sizes;
};

#endif  // CLOTHO_STATUS_BUFFER_HPP_
