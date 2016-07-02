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
#ifndef CLOTHO_CROSSOVER_GENERATOR_HPP_
#define CLOTHO_CROSSOVER_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/utility/bit_helper.hpp"
#include "clotho/utility/debruijn_bit_walker.hpp"
#include "clotho/data_spaces/crossover/crossover_method.hpp"
#include "clotho/data_spaces/generators/crossover_event_generator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class Crossover {
public:
    typedef GeneticSpaceType    genetic_space_type;

    typedef typename genetic_space_type::allele_type        allele_type;
    typedef typename genetic_space_type::association_type   association_type;

    typedef crossover_method< RNG, allele_type, association_type >  method_type;

//    typedef typename genetic_space_type::block_type         block_type;
//    typedef typename allele_type::position_type             position_type;

//    typedef typename genetic_space_type::sequence_iterator      sequence_iterator;
//    typedef typename genetic_space_type::genome_iterator        genome_iterator;
    typedef typename method_type::sequence_iterator             sequence_iterator;
    typedef typename method_type::genome_iterator               genome_iterator;
    typedef typename genetic_space_type::individual_id_type     individual_id_type;

//    typedef crossover_event_generator< RNG, position_type >     event_generator_type;
//
//    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
//    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    typedef std::vector< std::pair< individual_id_type, individual_id_type > >  mate_pair_type;
    typedef typename mate_pair_type::iterator                                   iterator;

    Crossover( RNG * rng, boost::property_tree::ptree & config ) :
        m_method( rng, config )
    { }

    void update( genetic_space_type * parental_genomes, mate_pair_type & parents, genetic_space_type * offspring_genomes ) {
        iterator mate_it = parents.begin(), mate_end = parents.end();

        m_method.initAlleles( parental_genomes->getAlleleSpace() );

//        m_event_gen.update( parental_genomes->getAlleleSpace().position_begin(), parental_genomes->getAlleleSpace().position_end() );
//
        offspring_genomes->getSequenceSpace().clear();

        size_t i = 0;
        size_t p_step = parental_genomes->getSequenceSpace().block_column_count();
        size_t c_step = offspring_genomes->getSequenceSpace().block_column_count();
        while( mate_it != mate_end ) {

            genome_iterator start = parental_genomes->begin_genome( mate_it->first );
            genome_iterator end = parental_genomes->end_genome( mate_it->first );
            sequence_iterator c = offspring_genomes->begin_sequence( i++ );

            m_method.crossover( start, end, c, p_step, c_step );
            
            start = parental_genomes->begin_genome( mate_it->second );
            end = parental_genomes->end_genome( mate_it->second );
            c = offspring_genomes->begin_sequence( i++ );

            m_method.crossover( start, end, c, p_step, c_step );

            ++mate_it;
        }
    }

    virtual ~Crossover() {}

protected:
    method_type m_method;
};

}   // namespace genetics
}   // namespace clotho

#include "clotho/data_spaces/crossover/crossover_spec.hpp"

#endif  // CLOTHO_CROSSOVER_GENERATOR_HPP_
