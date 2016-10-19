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
#ifndef CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_HPP_
#define CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_HPP_

#include "clotho/data_spaces/task/task.hpp"
#include "clotho/data_spaces/free_space/free_space_task.hpp"
#include "clotho/data_spaces/population_space/population_space.hpp"

#include <set>

namespace clotho {
namespace genetics {

template < class BlockType, class WeightType >
class free_space_task < population_space< BlockType, WeightType > > : public task {
public:
    typedef free_space_task< population_space< BlockType, WeightType > > self_type;

    typedef population_space< BlockType, WeightType >       space_type;
    typedef typename space_type::base_genome_type::sequence_type::block_type block_type;

    typedef typename space_type::genome_type                genome_type;

//    typedef typename space_type::base_genome_type *         genome_type;
    typedef typename space_type::genome_iterator            iterator;

    free_space_task( space_type * ss, block_type * fixed_dest, block_type * var_dest, unsigned int block_start, unsigned int block_end ) :
        m_space(ss)
        , m_destF( fixed_dest )
        , m_destV( var_dest )
        , m_start( block_start )
        , m_end( block_end )
    {}

    free_space_task( const self_type & other ) :
        m_space( other.m_space )
        , m_destF( other.m_destF )
        , m_destV( other.m_destV )
        , m_start( other.m_start )
        , m_end( other.m_end )
    {}

//    void operator()() {
//        iterator first = m_space->begin_genomes(), last = m_space->end_genomes();
//        unsigned int N = 0;
//
//        unsigned int i = 0;
//        while( first != last ) {
//            process( first->first );
//            N += first->second;
//
//            ++i;
//            ++first;
//        }
//
//#ifdef DEBUGGING
//        std::cerr << "Scanned " << i << " genomes" << std::endl;
//        BOOST_LOG_TRIVIAL(debug) << "Processed " << processed.size() << " sequences from " << N << " individuals";
//#endif  // DEBUGGING
//
//        assert( N == m_space->haploid_genome_count() );
//    }

//    void operator()() {
//        individual_iterator first = m_space->begin_individual(), last = m_space->end_individual();
//
//        std::set< base_genome_type * > processed;
//
//        while( first != last ) {
//            if( processed.find( first->first.get() ) == processed.end() ) {        
//                process( first->first );
//                processed.insert( first->first.get() );
//            }
//
//            if( processed.find( first->second.get() ) == processed.end() ) {
//                process( first->second );
//               processed.insert( first->second.get() );
//            }
//
//            ++first;
//        }
//    }
//
    void operator()() {
        typename space_type::individual_iterator first = m_space->begin_individual(), last = m_space->end_individual();
        unsigned int N = 0;
        while( first != last ) {
            process( first->first );
            process( first->second );
            ++first;
            N += 2;
        }

        std::cerr << "Scanned " << N << " haploid genomes" << std::endl;
    }

    void process ( genome_type hap ) {
        typedef typename space_type::base_genome_type::sequence_type::const_sequence_iterator seq_iterator;

        seq_iterator first, last;

        if( !hap ) {
            first = last;
        } else if( m_start >= hap->sequence_block_length() ) {
            first = hap->end_sequence();
            last = hap->end_sequence();
        } else {
            first = hap->begin_sequence() + m_start;
            last = (( m_end >= hap->sequence_block_length() ) ? hap->end_sequence() : hap->begin_sequence() + m_end);
        }

        block_type * df = m_destF + m_start, * de = m_destF + m_end;
        block_type * dv = m_destV + m_start;

        while( first != last ) {
            assert( df != de );

            const block_type b = *first++;
            *df++ &= b;
            *dv++ |= b;
        }

//        static const block_type ALL_UNSET = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
//
        //std::cerr << "zeroing fixed length: " << (de - df) << std::endl;
        while( df != de ) {
            *df++ = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
        }
    }

//    void process( const genome_type & hap ) {
//        typedef typename space_type::base_genome_type::sequence_type::const_sequence_iterator seq_iterator;
//
//        block_type * df = m_destF + m_start, * de = m_destF + m_end;
//        if( !hap || m_start >= hap->sequence_block_length() ) {
////            unsigned int i = m_start;
////            while( i < m_end ) {
////                m_destF[ i++ ] = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
////            }
//            while( df != de ) {
//                *df++ = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
//            }
//        } else {
//            seq_iterator first = hap->begin_sequence() + m_start, 
//                last = (( m_end >= hap->sequence_block_length() ) ? hap->end_sequence() : hap->begin_sequence() + m_end);
//
////            unsigned int i = m_start;
////            while( first != last ) {
////                const block_type b = *first++;
////                m_destF[ i ] &= b;
////                m_destV[ i++ ] |= b;
////            }
////
////            while( i < m_end ) {
////                m_destF[ i++ ] = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
////            }
//
//            block_type * dv = m_destV + m_start;
//            while( first != last ) {
//                const block_type b = *first++;
//                *df++ &= b;
//                *dv++ |= b;
//            }
//
//            while( df != de ) {
//                *df++ = space_type::base_genome_type::sequence_type::bit_helper_type::ALL_UNSET;
//            }
//        }
//    }

    virtual ~free_space_task() {}

protected:
    space_type  * m_space;

    block_type  * m_destF, * m_destV;
    unsigned int m_start, m_end;
};

}   // namespace genetics
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_TASK_POPULATION_SPACE_HPP_
