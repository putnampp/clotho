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
#ifndef DEVICE_SEQUENCE_SPACE_HELPER_HPP_
#define DEVICE_SEQUENCE_SPACE_HELPER_HPP_

template < class IntType >
void dump_data( std::ostream & out, const device_sequence_space< IntType > & space ) {
    if( space.seq_count == 0 || space.seq_width == 0 ) {
        out << "[]";
        return;
    }

    typedef typename device_sequence_space< IntType >::int_type int_type;

    int_type * dSeq;
    int_type * seqs = new int_type[ space.size ];

    assert( cudaMalloc( (void **) &dSeq, space.size * sizeof( int_type ) ) == cudaSuccess );
    copy_heap<<< 1, 1 >>>( dSeq, space.sequences, space.size );
    cudaDeviceSynchronize();

    assert( cudaMemcpy( seqs, dSeq, space.size * sizeof( int_type ), cudaMemcpyDeviceToHost ) == cudaSuccess );

    int_type * ptr = seqs;
    out << std::hex;
    out << "[";
    for( int i = 0; i < space.seq_count; ++i ) {
        if( i ) {
            out << ",\n";
        }
        out << "[0x" << std::setw( sizeof(int_type) * 2) << std::setfill('0') << (*ptr++);
        for( int j = 1; j < space.seq_width; ++j ) {
            out << ",0x" << std::setw( sizeof(int_type) * 2) << std::setfill('0') << (*ptr++);
        }
        out << "]";
    }

    out << "]";
    out << std::dec;

    cudaFree( dSeq );
    delete seqs;
}

template < class IntType >
std::ostream & operator<<( std::ostream & out, const device_sequence_space< IntType > & space ) {
    out << "{";
    out << "\n\"data\": ";
#ifdef DUMP_HEAP_VARIABLES 
    dump_data( out, space );
#else
    out << "[]";
#endif  // DUMP_HEAP_VARIABLES

    out << ",\n\"seq_count\": " << space.seq_count;
    out << ",\n\"seq_width\": " << space.seq_width;
    out << ",\n\"size\": " << space.size;
    out << ",\n\"capacity\": " << space.capacity;
    out << "}";

    return out;
}

namespace clotho {
namespace utility {

template < >
void get_state( boost::property_tree::ptree & state, const device_sequence_space< unsigned int > & obj ) {
    state.put( "dimensions.rows", obj.seq_count );
    state.put( "dimensions.columns", obj.seq_width );

    state.put( "size", obj.size);
    state.put( "capacity", obj.capacity );
}

}   // namespace utility
}   // namespace clotho

#endif  // DEVICE_SEQUENCE_SPACE_HELPER_HPP_
