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
#ifndef CLOTHO_FREE_SPACE_HPP_
#define CLOTHO_FREE_SPACE_HPP_

//#include "clotho/utility/bit_helper.hpp"
//#include "clotho/utility/debruijn_bit_walker.hpp"

#include "clotho/utility/state_object.hpp"
#include "clotho/utility/log_helper.hpp"

#include "clotho/data_spaces/free_space/free_space_evaluator.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
class FreeSpaceAnalyzer {
public:
    typedef GeneticSpaceType    genetic_space_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef clotho::genetics::free_space_evaluator< association_type > evaluator_type;

//    typedef typename association_type::raw_block_pointer    block_pointer;
//
//    typedef typename association_type::block_type           block_type;

    typedef size_t *                                    index_vector;
    typedef size_t *                                    iterator;
    typedef size_t *                                    const_iterator;

//    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
//
//    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    FreeSpaceAnalyzer():
        m_indices( NULL )
        , m_fixed_count(0)
        , m_lost_count(0)
        , m_free_count(0)
        , m_width(0)
        , m_size(0)
    {}

    FreeSpaceAnalyzer( genetic_space_type & gs ) :
        m_indices( NULL )
        , m_fixed_count(0)
        , m_lost_count(0)
        , m_free_count(0)
        , m_width(0)
        , m_size(0)
    {
        update( gs );
    }

    void update( genetic_space_type & gs ) {
        size_t M = gs.allele_count();
        resize( M ) ;

        m_free_count = 0;
        m_fixed_count = m_width;
        m_lost_count = 2 * m_width;

//        evaluator_type eval;
        eval( gs.getSequenceSpace(), m_indices, m_fixed_count, m_lost_count, m_free_count, m_width );

        m_fixed_count -= m_width;
        m_lost_count -= (2 * m_width);

        assert( m_fixed_count + m_lost_count == m_free_count );
    }

    size_t free_size() const {
        return m_free_count;
    }

    iterator free_begin() {
        return m_indices;
    }

    iterator free_end() {
        return m_indices + m_free_count;
    }

    const_iterator free_begin() const {
        return m_indices;
    }

    const_iterator free_end() const {
        return m_indices + m_free_count;
    }

    size_t fixed_size() const {
        return m_fixed_count;
    }

    iterator fixed_begin() {
        return m_indices + m_width;
    }

    iterator fixed_end() {
        return m_indices + m_width + m_fixed_count;
    }

    const_iterator fixed_begin() const {
        return m_indices + m_width;
    }

    const_iterator fixed_end() const {
        return m_indices + m_width + m_fixed_count;
    }

    size_t lost_size() const {
        return m_lost_count;
    }

    iterator lost_begin() {
        return m_indices + 2 * m_width;
    }

    iterator lost_end() {
        return m_indices + 2 * m_width + m_lost_count;
    }

    const_iterator lost_begin() const {
        return m_indices + 2 * m_width;
    }

    const_iterator lost_end() const {
        return m_indices + 2 * m_width + m_lost_count;
    }

    virtual ~FreeSpaceAnalyzer() {
        if( m_indices != NULL ) {
            delete [] m_indices;
        }
    }

protected:
    
    void resize( size_t s ) {
        if( s > m_width ) {
            if( m_indices != NULL ) {
                delete [] m_indices;
            }

            m_indices = new size_t[ 3 * s ];
            m_size = 3 * s;
            m_width = s;
        }
    }

    size_t  * m_indices;
    size_t  m_fixed_count, m_lost_count, m_free_count;
    size_t  m_width, m_size;
    evaluator_type eval;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class GeneticSpaceType >
struct state_getter< clotho::genetics::FreeSpaceAnalyzer< GeneticSpaceType > > {
    typedef clotho::genetics::FreeSpaceAnalyzer< GeneticSpaceType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        boost::property_tree::ptree fr;
        clotho::utility::add_value_array( fr, obj.free_begin(), obj.free_end() );
        
        s.put_child( "free", fr );
    }
};

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_HPP_
