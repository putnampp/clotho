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
#include "clotho/data_spaces/generators/crossover_event_generator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class GeneticSpaceType >
class Crossover {
public:
    typedef GeneticSpaceType    genetic_space_type;

    typedef typename genetic_space_type::allele_type        allele_type;
    typedef typename genetic_space_type::block_type         block_type;
    typedef typename allele_type::position_type             position_type;

    typedef typename genetic_space_type::sequence_iterator  sequence_iterator;
    typedef typename genetic_space_type::genome_iterator    genome_iterator;

    typedef crossover_event_generator< RNG, position_type >       event_generator_type;

    typedef clotho::utility::debruijn_bit_walker< block_type >  bit_walker_type;
    typedef clotho::utility::BitHelper< block_type >            bit_helper_type;

    typedef std::vector< std::pair< size_t, size_t > >  mate_pair_type;
    typedef typename mate_pair_type::iterator           iterator;

    Crossover( RNG * rng, boost::property_tree::ptree & config ) :
        m_event_gen(rng, config )
    {}

    void update( genetic_space_type * parental_genomes, mate_pair_type & parents, genetic_space_type * offspring_genomes ) {
        iterator mate_it = parents.begin(), mate_end = parents.end();

        m_event_gen.update( parental_genomes->getAlleleSpace().position_begin(), parental_genomes->getAlleleSpace().position_end() );

        size_t i = 0;
        while( mate_it != mate_end ) {
            genome_iterator p = parental_genomes->getGenomeAt( mate_it->first );
            sequence_iterator c = offspring_genomes->getSequenceAt( i++ );

            crossover( p, c, parental_genomes->getAlleleSpace() );
            
            genome_iterator q = parental_genomes->getGenomeAt( mate_it->second );
            c = offspring_genomes->getSequenceAt( i++ );
            crossover( q, c, parental_genomes->getAlleleSpace() );

            ++mate_it;
        }
    }

    virtual ~Crossover() {}

protected:
    void crossover( genome_iterator & p, sequence_iterator & s, allele_type & all ) {
        m_event_gen.generate();
        size_t j = 0;
        while( p.hasNext() ) {
            std::pair< block_type, block_type > g = p.next_pair();

            block_type hets = g.first ^ g.second;
            block_type sec  = bit_helper_type::ALL_UNSET;   // mask state from second strand

            while( hets ) {
                size_t idx = bit_walker_type::unset_next_index( hets ) + j;

                if( m_event_gen( idx ) ) {
                    sec |= bit_helper_type::bit_offset( idx );
                }
            }

            hets = ((g.first & ~sec) | (g.second & sec));
            s.write_next( hets );

            j += bit_helper_type::BITS_PER_BLOCK;
        }
    }

    event_generator_type    m_event_gen;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_CROSSOVER_GENERATOR_HPP_
