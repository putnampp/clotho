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

#include "clotho/data_spaces/free_space/free_space_details.hpp"
#include "clotho/data_spaces/free_space/free_space_evaluator.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType, class SizeType = unsigned int >
class FreeSpaceAnalyzer : public free_space_details< SizeType > {
public:
    typedef free_space_details< SizeType >                  base_type;
    typedef GeneticSpaceType    genetic_space_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef typename base_type::size_type                   size_type;
    typedef clotho::genetics::free_space_evaluator< association_type, size_type > evaluator_type;

//    typedef typename association_type::raw_block_pointer    block_pointer;
//
//    typedef typename association_type::block_type           block_type;

    typedef typename base_type::index_vector                index_vector;
    typedef typename base_type::iterator                    iterator;
    typedef typename base_type::const_iterator              const_iterator;


//    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
//
//    typedef clotho::utility::BitHelper< block_type >    bit_helper_type;

    FreeSpaceAnalyzer() {}

    //FreeSpaceAnalyzer( genetic_space_type & gs ) {
    //    update( gs );
    //}

    void update( genetic_space_type & gs ) {
        size_type M = gs.allele_count();
        this->resize( M ) ;

        this->m_free_count = 0;
        this->m_fixed_count = this->m_width;
        this->m_lost_count = 2 * this->m_width;

        eval( gs.getSequenceSpace(), this->m_indices, this->m_fixed_count, this->m_lost_count, this->m_free_count, this->m_width );

        this->m_fixed_count -= this->m_width;
        this->m_lost_count -= (2 * this->m_width);

        assert( this->m_fixed_count + this->m_lost_count == this->m_free_count );
    }

    virtual ~FreeSpaceAnalyzer() {}

protected:
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
        s.add( "free_size", obj.free_size() );
        s.add( "variable_count", obj.variable_count() );
    }
};

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_FREE_SPACE_HPP_
