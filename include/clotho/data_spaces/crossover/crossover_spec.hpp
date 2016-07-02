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
#ifndef CLOTHO_CROSSOVER_GENERATOR_SPEC_HPP_
#define CLOTHO_CROSSOVER_GENERATOR_SPEC_HPP_

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

    typedef typename method_type::row_vector                    row_vector;
    
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

        while( mate_it != mate_end ) {

            row_vector a = parental_genomes->getSequenceSpace().getRow( 2 * mate_it->first );
            row_vector b = parental_genomes->getSequenceSpace().getRow( 2 * mate_it->first + 1 );
            row_vector c = m_method.crossover( a, b );

            offspring_genomes->getSequenceSpace().push_back( c );
            
            row_vector d = parental_genomes->getSequenceSpace().getRow( 2 * mate_it->second );
            row_vector e = parental_genomes->getSequenceSpace().getRow( 2 * mate_it->second + 1 );
            row_vector f = m_method.crossover( d, e );

            offspring_genomes->getSequenceSpace().push_back( f );

            ++mate_it;
        }

//        offspring_genomes->getSequenceSpace().finalize();
    }

    virtual ~Crossover() {}

protected:
    method_type m_method;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_GENERATOR_SPEC_HPP_
