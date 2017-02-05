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
#ifndef CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_
#define CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_

#include "clotho/data_spaces/free_space/free_space_buffer.hpp"
#include "clotho/data_spaces/free_space/free_space_details.hpp"

namespace clotho {
namespace genetics {

template < class BlockType, class SizeType >
class free_space_accumulator : public free_space_details< SizeType >, public free_space_buffer< BlockType > {
public:
    typedef free_space_accumulator< BlockType, SizeType >   self_type;
    typedef free_space_details< SizeType >                  base_type;
    typedef free_space_buffer< BlockType >                  buffer_type;

    free_space_accumulator( ) {}

    void resetBuffers( buffer_type * other ) {
        this->reset( other->size() );
        updateBuffers( other );
    }

    void updateBuffers( buffer_type * buf ) {
        this->updateVariableBuffer( buf->getVariableBuffer(), buf->size());
        this->updateFixedBuffer( buf->getFixedBuffer(), buf->size() );
    }

    void finalize( unsigned int max_alleles ) {
        this->resize( max_alleles );
        this->analyze_free_indices( this->getFixedBuffer(), this->getVariableBuffer(), this->size(), max_alleles );
    }

    virtual ~free_space_accumulator() { }

protected:
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class BlockType, class SizeType >
struct state_getter< clotho::genetics::free_space_accumulator< BlockType, SizeType > > {
    typedef clotho::genetics::free_space_accumulator< BlockType, SizeType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree fr;
        clotho::utility::add_value_array( fr, obj.free_begin(), obj.free_end() );

        s.put_child( "free", fr );
        s.add( "free_size", obj.free_size() );
        s.add( "variable_size", obj.variable_count() );
    }
};

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_ACCUMULATOR_HPP_
