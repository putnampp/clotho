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
#ifndef CLOTHO_ROW_CROSSOVER_HPP_
#define CLOTHO_ROW_CROSSOVER_HPP_


#include <boost/property_tree/ptree.hpp>

#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"
#include "clotho/data_spaces/crossover/crossover_method.hpp"
#include "clotho/data_spaces/generators/crossover_event_generator.hpp"

#include "clotho/data_spaces/crossover/crossover.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/association_matrix/row_grouped_association_matrix.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class AlleleSpaceType, class BlockType >
class Crossover< RNG, genetic_space< AlleleSpaceType, BlockType, row_grouped< 1 > > > {
public:
    typedef genetic_space< AlleleSpaceType, BlockType, row_grouped< 1 > > genetic_space_type;

    typedef typename genetic_space_type::allele_type        allele_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef crossover_method< RNG, allele_type, association_type >  method_type;

    typedef typename method_type::sequence_iterator             sequence_iterator;
    
    typedef typename genetic_space_type::individual_id_type     individual_id_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
    typedef typename mate_pair_type::iterator                                   iterator;

    Crossover( RNG * rng, boost::property_tree::ptree & config ) :
        m_method( rng, config )
    {
#ifdef DEBUGGING
        std::cerr << "Using specialized crossover method" << std::endl;
#endif  // DEBUGGING
    }

    void update( genetic_space_type * parental_genomes, mate_pair_type & parents, genetic_space_type * offspring_genomes ) {
        iterator mate_it = parents.begin(), mate_end = parents.end();

        m_method.initAlleles( parental_genomes->getAlleleSpace() );

        size_t i = 0;
        while( mate_it != mate_end ) {

            size_t idx = 2 * mate_it->first;
            sequence_iterator first = parental_genomes->begin_sequence( idx );
            unsigned int N = parental_genomes->getSequenceSpace().getSoftSize( idx++ );

            sequence_iterator second = parental_genomes->begin_sequence( idx );
            unsigned int M = parental_genomes->getSequenceSpace().getSoftSize( idx );

            sequence_iterator c = offspring_genomes->begin_sequence( i );
            unsigned int W = m_method.crossover( first, N, second, M, c );
            offspring_genomes->getSequenceSpace().updateSoftSize( i, W );

#ifdef  DEBUGGING
            assert( N <= parental_genomes->getSequenceSpace().block_column_count());
            assert( M <= parental_genomes->getSequenceSpace().block_column_count());
            assert( W <= offspring_genomes->getSequenceSpace().block_column_count());
#endif  // DEBUGGING
            ++i;

            idx = 2 * mate_it->second;
            first = parental_genomes->begin_sequence( idx );
            N = parental_genomes->getSequenceSpace().getSoftSize( idx++ );

            second = parental_genomes->begin_sequence( idx );
            M = parental_genomes->getSequenceSpace().getSoftSize( idx );

            c = offspring_genomes->begin_sequence( i );
            W = m_method.crossover( first, N, second, M, c );
            offspring_genomes->getSequenceSpace().updateSoftSize(i, W );

#ifdef  DEBUGGING
            assert( N <= parental_genomes->getSequenceSpace().block_column_count());
            assert( M <= parental_genomes->getSequenceSpace().block_column_count());
            assert( W <= offspring_genomes->getSequenceSpace().block_column_count());
#endif  // DEBUGGING
            ++i;

            ++mate_it;
        }
    }

    virtual ~Crossover() {}

protected:
    method_type m_method;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ROW_CROSSOVER_HPP_
