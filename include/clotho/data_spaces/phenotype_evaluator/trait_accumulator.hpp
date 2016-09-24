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
#ifndef CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_
#define CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_

#include <vector>

#include "clotho/data_spaces/growable2D.hpp"
#include "clotho/data_spaces/trait_helper_of.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"

#include "clotho/utility/state_object.hpp"

#include "clotho/data_spaces/phenotype_evaluator/weight_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class GeneticSpaceType >
class TraitWeightAccumulator : public growable2D {
public:
    typedef TraitWeightAccumulator< GeneticSpaceType >          self_type;

    typedef GeneticSpaceType                                                    genetic_space_type;
    typedef typename genetic_space_type::block_type                             block_type;
    typedef typename genetic_space_type::association_type                       association_type;
    typedef typename association_type::raw_block_pointer    block_pointer;

    typedef typename genetic_space_type::allele_type            allele_type;
    typedef typename allele_type::weight_type                   weight_type;
    typedef typename allele_type::weight_pointer                weight_pointer;
    
    typedef clotho::genetics::trait_helper_of< allele_type >    trait_helper_type;
    typedef typename trait_helper_type::type                    trait_vector_type;

    typedef std::vector< trait_vector_type >                    accumulator_type;

    typedef clotho::genetics::weight_accumulator< allele_type, association_type > evaluator_type;

    TraitWeightAccumulator( ) :
        m_trait_weights(NULL)
        , m_rows(0)
        , m_trait_count(0)
        , m_size(0)
    {}

    TraitWeightAccumulator( genetic_space_type & genomes ) :
        m_trait_weights(NULL)
        , m_rows(0)
        , m_trait_count(0)
        , m_size(0)
    {
        update( genomes );
    }

    void update( genetic_space_type & genomes ) {

        this->grow( genomes.sequence_count(), genomes.getAlleleSpace().trait_count() );

        memset( m_trait_weights, 0, m_size * sizeof(weight_type) );

        if( genomes.getAlleleSpace().isAllNeutral() ) {
            return;
#ifdef DEBUGGING
        } else {
            BOOST_LOG_TRIVIAL(debug) << "Not all alleles are neutral!";
#endif // DEBUGGING
        }

        evaluator_type eval;

        eval( genomes.getAlleleSpace(), genomes.getSequenceSpace(), m_trait_weights, m_trait_count );
    }

//    std::shared_ptr< trait_vector_type >   getTraitAt( size_t idx ) {
//        assert( 0 <= idx && idx < m_rows );
//
//       idx *= m_trait_count;
//
//        std::shared_ptr< trait_vector_type>  t( new trait_vector_type( m_trait_weights + idx, m_trait_weights + idx + m_trait_count ));
//
//       return t;
//    }

    weight_pointer begin_trait_weight( size_t idx ) {
        assert( 0 <= idx && idx < m_rows );
        return m_trait_weights + idx * m_trait_count;
    }

    weight_pointer end_trait_weight( size_t idx ) {
        assert( 0 <= idx && idx < m_rows );
        return m_trait_weights + (idx + 1) * m_trait_count;
    }

    size_t size() const {
        return m_rows;
    }

    size_t trait_count() const {
        return m_trait_count;
    }

    virtual size_t grow( size_t seq ) {
        size_t t = trait_count();
        if( t == 0 ) {
            t = 1;
        }

        this->resize( seq, t );
        return m_rows;
    }

    virtual size_t grow( size_t seq, size_t traits ) {
        this->resize( seq, traits );
        return m_rows;
    }

    virtual ~TraitWeightAccumulator() {
        if( m_trait_weights != NULL ) {
            delete [] m_trait_weights;
        }
    }

protected:

    virtual void resize( size_t rows, size_t columns ) {
        size_t new_size = rows * columns;

        assert( 0 <= new_size );

        if( m_size < new_size ) {
            if( m_trait_weights != NULL ) {
                delete [] m_trait_weights;
            }
            m_trait_weights = new weight_type[ new_size ];
            m_size = new_size;
        }

        m_rows = rows;
        m_trait_count = columns;
    }

    weight_pointer      m_trait_weights;
    size_t              m_rows, m_trait_count,m_size;
};

}   // namespace genetics
}   // namespace clotho

namespace clotho {
namespace utility {

template < class GeneticSpaceType >
struct state_getter< clotho::genetics::TraitWeightAccumulator< GeneticSpaceType > > {
    typedef clotho::genetics::TraitWeightAccumulator< GeneticSpaceType > object_type;

    void operator()( boost::property_tree::ptree & s, object_type & obj ) {
        size_t i = 0, N = obj.size();
        while( i < N ) {
//            std::shared_ptr< typename object_type::trait_vector_type > t = obj.getTraitAt(i++);
            boost::property_tree::ptree p;
            
            clotho::utility::add_value_array( p, obj.begin_trait_weight(i), obj.end_trait_weight(i) );

            s.push_back( std::make_pair( "", p ) );
            ++i;
        }
    
    }
};


}   // namespace utility
}   // namespace clotho

#endif  // CLOTHO_TRAIT_WEIGHT_ACCUMULATOR_HPP_
