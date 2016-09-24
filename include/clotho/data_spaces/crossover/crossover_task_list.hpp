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
#ifndef CLOTHO_CROSSOVER_TASK_LIST_HPP_
#define CLOTHO_CROSSOVER_TASK_LIST_HPP_

#include "clotho/data_spaces/crossover/tail_crossover_task.hpp"
#include "clotho/data_spaces/crossover/segment_crossover_task.hpp"
#include "clotho/data_spaces/crossover/copy_crossover_task.hpp"
#include "clotho/data_spaces/crossover/block_crossover_and_mutate_task.hpp"
#include "clotho/data_spaces/mutation/block_mutate_task.hpp"

#include <boost/random/bernoulli_distribution.hpp>

namespace clotho {
namespace genetics {

template < class Classifier, class BlockType >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res ) {
    typedef copy_crossover_task< BlockType > copy_type;
    typedef segment_crossover_task< Classifier, BlockType > segment_type;
    typedef tail_crossover_task< Classifier, BlockType, top_strand_tail > top_tail_type;
    typedef tail_crossover_task< Classifier, BlockType, bottom_strand_tail > bottom_tail_type;

    if( s0_len == 0 && s1_len == 0 ) return;

    if( cls.event_count() == 0 ) {
        // there are no crossover cls
        // therefore, offspring strand will be a copy of the top strand
        copy_type t( s0, res, s0_len );
        t();
    } else if( s0_len < s1_len ) {
        // top strand is shorter than bottom strand
        segment_type s( cls, s0, s1, res, 0, s0_len );
        s();

        bottom_tail_type t( cls, s1, res, s0_len, s1_len );
        t();
    } else if( s1_len < s0_len ) {
        // bottom strand is shorter than top strand
        segment_type s( cls, s0, s1, res, 0, s1_len );
        s();

        top_tail_type t( cls, s0, res, s1_len, s0_len );
        t();
    } else {
        // both strands are equal in length
        segment_type s( cls, s0, s1, res, 0, s0_len );
        s();
    }
}

template < class Classifier, class BlockType, class PoolType >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service ) {
    typedef copy_crossover_task< BlockType > copy_type;
    typedef segment_crossover_task< Classifier, BlockType > segment_type;
    typedef tail_crossover_task< Classifier, BlockType, top_strand_tail > top_tail_type;
    typedef tail_crossover_task< Classifier, BlockType, bottom_strand_tail > bottom_tail_type;

    if( s0_len == 0 && s1_len == 0 ) return;

    if( cls.event_count() == 0 ) {
        // there are no crossover cls
        // therefore, offspring strand will be a copy of the top strand
        copy_type c( s0, res, s0_len );
        service.post( c );
    } else if( s0_len < s1_len ) {
        // top strand is shorter than bottom strand
        segment_type s( cls, s0, s1, res, 0, s0_len );
        service.post( s );

        bottom_tail_type b( cls, s1, res, s0_len, s1_len );
        service.post( b );
    } else if( s1_len < s0_len ) {
        // bottom strand is shorter than top strand
        segment_type s( cls, s0, s1, res, 0, s1_len );
        service.post( s );

        top_tail_type t( cls, s0, res, s1_len, s0_len );
        service.post( t );
    } else {
        // both strands are equal in length
        segment_type s( cls, s0, s1, res, 0, s0_len );
        service.post( s );
    }
}

template < class Classifier, class BlockType, class PoolType >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service, unsigned int start_idx, unsigned int end_idx ) {
    typedef copy_crossover_task< BlockType > copy_type;
    typedef segment_crossover_task< Classifier, BlockType > segment_type;
    typedef tail_crossover_task< Classifier, BlockType, top_strand_tail > top_tail_type;
    typedef tail_crossover_task< Classifier, BlockType, bottom_strand_tail > bottom_tail_type;

    if( s0_len == 0 && s1_len == 0 ) return;

    if( start_idx >= end_idx ) return;

    // start_idx < end_idx

    if( cls.event_count() == 0 ) {
        // there are no crossover cls
        // therefore, offspring strand will be a copy of the top strand
        unsigned int len = ((end_idx < s0_len) ? end_idx : s0_len);

        // catches whether s0_len <= start_idx
        if( len <= start_idx ) return;

        len -= start_idx;

        copy_type c( s0 + start_idx, res + start_idx, len );
        service.post( c );
    } else if( start_idx < s0_len ) {
        if( start_idx < s1_len ) {
            // start_idx < s0_len && start_idx < s1_len
            if( end_idx <= s0_len) {
                if ( end_idx <= s1_len ) {
                    // end_idx <= s0_len && end_idx <= s1_len
                    // end_idx - start_idx > 0 results from start_idx < end_idx
                    segment_type s( cls, s0, s1, res, start_idx, (end_idx - start_idx) );
                    service.post( s );
                } else {
                    // s1_len < end_idx <= s0_len
                    // s1_len - start_idx > 0 results from start_idx < s1_len
                    // end_idx - s1_len > 0 results from s1_len < end_idx
                    segment_type s( cls, s0, s1, res, start_idx, (s1_len - start_idx) );
                    service.post( s );

                    top_tail_type t( cls, s0, res, s1_len, (end_idx - s1_len));
                    service.post( t );
                }
            } else if( end_idx <= s1_len ) {
                // s0_len < end_idx <= s1_len
                // s0_len - start_idx > 0 results from start_idx < s0_len
                // end_idx - s0_len > 0 results from !(end_idx <= s0_len)
                segment_type s( cls, s0, s1, res, start_idx, (s0_len - start_idx) );
                service.post( s );

                bottom_tail_type b( cls, s1, res, s0_len, (end_idx - s0_len));
                service.post( b );
            } else if( s1_len < s0_len ) {
                // end_idx > s0_len && end_idx > s1_len && s0_len < s1_len
                // s0_len - start_idx > 0 results from (start_idx < s0_len)
                // s1_len - s0_len > 0 results from (s1_len < s0_len)
                segment_type s( cls, s0, s1, res, start_idx, (s0_len - start_idx) );
                service.post( s );

                bottom_tail_type b( cls, s1, res, s0_len, ( s1_len - s0_len ) );
                service.post( b );
            } else if( s0_len < s1_len ) {
                // end_idx > s0_len && end_idx > s1_len && s1_len < s0_len
                // s1_len - start_idx > 0 results from (start_idx < s1_len)
                // s0_len - s1_len > 0 results from (s1_len >= s0_len) && (s0_len < s1_len)
                segment_type s( cls, s0, s1, res, start_idx, (s1_len - start_idx) );
                service.post( s );

                top_tail_type t( cls, s1, res, s1_len, ( s0_len - s1_len ) );
                service.post( t );
            } else {
                // end_idx > s0_len && end_idx > s1_len && s1_len == s0_len
                // s1_len - start_idx > 0 results from (start_idx < s1_len)
                segment_type s( cls, s0, s1, res, start_idx, (s1_len - start_idx) );
                service.post( s );
            }
        } else  {
            // s1_len <= start_idx < s0_len
            unsigned int e = ((end_idx < s0_len) ? end_idx : s0_len );
            top_tail_type t( cls, s0, res, start_idx, (e - start_idx) );
            service.post( t ); 
        }
    } else if( start_idx < s1_len ) {
        // s0_len <= start_idx < s1_len
        unsigned int e = ((end_idx < s1_len ) ? end_idx : s1_len );

        // e -start_idx > 0 results from (start_idx < end_idx) && (start_idx < s1_len)
        bottom_tail_type b( cls, s1, res, start_idx, (e - start_idx) );
        service.post( b );
    }
}

template < class Classifier, class BlockType, class PoolType, class RNG >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service, double strand_bias, RNG & rng ) {
    boost::random::bernoulli_distribution< double > dist( strand_bias );

    if( dist(rand) ) {
        // consider s1 as top strand, s0 as bottom strand
        make_crossover_tasks( cls, s1, s1_len, s0, s0_len, res, service );
    } else {
        make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service );
    }
}

template < class Classifier, class BlockType, class PoolType >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service, bool should_swap_strands ) {

    if( should_swap_strands ) {
        // consider s1 as top strand, s0 as bottom strand
        make_crossover_tasks( cls, s1, s1_len, s0, s0_len, res, service );
    } else {
        make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service );
    }
}

template < class Classifier, class BlockType, class RNG >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, double strand_bias, RNG & rng ) {
    boost::random::bernoulli_distribution< double > dist( strand_bias );

    if( dist(rand) ) {
        // consider s1 as top strand, s0 as bottom strand
        make_crossover_tasks( cls, s1, s1_len, s0, s0_len, res );
    } else {
        make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res );
    }
}

template < class Classifier, class BlockType >
static void make_crossover_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, bool should_swap_strands ) {

    if( should_swap_strands ) {
        make_crossover_tasks( cls, s1, s1_len, s0, s0_len, res );
    } else {
        make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res );
    }
}
template < class Classifier, class BlockType, class PoolType >
static void make_crossover_and_mutate_task( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service, const block_mutate< BlockType > & mut, unsigned int idx ) {
    typedef clotho::utility::BitHelper< BlockType > bit_helper_type;
    typedef block_crossover_and_mutate_task< Classifier, BlockType > task_type;
    typedef block_mutate_task< BlockType > mutate_task_type;

    BlockType top = ((idx < s0_len) ? s0[ idx ] : bit_helper_type::ALL_UNSET );

    if( cls.event_count() == 0 ) {
        // no crossover events
        // therefore copy and mutate the top strand
        mutate_task_type m( top, res + idx, mut );
        service.post( m );
    } else {
        BlockType bottom = ((idx < s1_len ) ? s1[ idx ] : bit_helper_type::ALL_UNSET );
        task_type t( cls, top, bottom, res, idx, mut );

        service.post( t );
    }
}

template < class Classifier, class BlockType, class PoolType, class FreeIterator >
static void make_reproduction_tasks( const Classifier & cls, BlockType * s0, unsigned int s0_len, BlockType * s1, unsigned int s1_len, BlockType * res, PoolType & service, FreeIterator first, FreeIterator last ) {

    if( first == last ) {
        // no mutations
        // therefore follow standard crossover rules
        make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service );
    } else {
        
        typedef clotho::utility::BitHelper< BlockType > bit_helper_type;
        typedef block_mutate< BlockType >               mutate_type;

        FreeIterator tmp = first;
        unsigned int prev_block = *tmp++;
        prev_block /= bit_helper_type::BITS_PER_BLOCK;

        unsigned int start = 0;
        while( tmp != last ) {
            unsigned int cur_block = *tmp;
            cur_block /= bit_helper_type::BITS_PER_BLOCK;
            if( prev_block != cur_block ) {
                if( start != prev_block ) {
                    // schedule crossover from last mutation block to current mutation block
                    make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service, start, prev_block );
                }

                mutate_type m( first, tmp );
                make_crossover_and_mutate_task( cls, s0, s0_len, s1, s1_len, res, service, m, prev_block );

                first = tmp;

                start = prev_block + 1;
                prev_block = cur_block;
            }
            ++tmp;
        }

        if( first != tmp ) {
            if( start != prev_block ) {
                make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service, start, prev_block );
            }
            mutate_type m( first, tmp );
            make_crossover_and_mutate_task( cls, s0, s0_len, s1_len, res, service, m, prev_block );

            start = prev_block + 1;
        }

        prev_block = ((s0_len > s1_len) ? s0_len : s1_len );
        if( start < prev_block ) {
            make_crossover_tasks( cls, s0, s0_len, s1, s1_len, res, service, start, prev_block );
        }
    }
}

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_TASK_LIST_HPP_
