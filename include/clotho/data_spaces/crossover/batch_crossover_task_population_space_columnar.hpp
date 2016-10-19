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
#ifndef CLOTHO_BATCH_CROSSOVER_TASK_POPULATION_SPACE_COLUMNAR_HPP_
#define CLOTHO_BATCH_CROSSOVER_TASK_POPULATION_SPACE_COLUMNAR_HPP_

#include "clotho/data_spaces/population_space/population_space_columnar.hpp"
#include <boost/random/bernoulli_distribution.hpp>
#include "clotho/data_spaces/crossover/position_classifier.hpp"
#include "clotho/data_spaces/generators/position_distribution_helper.hpp"
#include "clotho/data_spaces/generators/crossover_event_distribution_helper.hpp"

#include "clotho/data_spaces/crossover/block_crossover.hpp"

#include "clotho/data_spaces/task/task.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class MatePairType, class BlockType, class WeightType, class AlleleSpaceType >
class batch_crossover_task< RNG, MatePairType, population_space_columnar< BlockType, WeightType > , AlleleSpaceType > : public task {
public:
    typedef batch_crossover_task< RNG, MatePairType, population_space_columnar< BlockType, WeightType >, AlleleSpaceType > self_type;

    typedef population_space_columnar< BlockType, WeightType >  space_type;
    typedef AlleleSpaceType                                     allele_type;
    typedef RNG                                                 random_engine_type;

    typedef typename space_type::block_type                     block_type;

    typedef typename space_type::individual_pointer             individual_pointer;
    typedef typename space_type::genome_pointer                 genome_pointer;


    typedef MatePairType                                        mate_pair_type;
    typedef typename mate_pair_type::iterator                   iterator;
    typedef typename mate_pair_type::const_iterator             const_iterator;

    typedef PositionClassifier< typename allele_type::position_vector > classifier_type;
    typedef typename classifier_type::event_type                        event_type;

    typedef block_crossover< classifier_type, BlockType >               crossover_type;
    typedef typename crossover_type::bit_helper_type                    bit_helper_type;

    typedef typename position_distribution_helper< typename allele_type::position_type >::type  position_distribution_type;
    typedef typename crossover_event_distribution_helper< double >::type                    event_distribution_type;

    batch_crossover_task( random_engine_type * rng, space_type * parents, space_type * offspring, allele_type * alleles, unsigned int off_idx, const_iterator first, const_iterator last, double recomb_rate, double seq_bias ) :
        m_rng( rng )
        , m_parent_pop( parents )
        , m_offspring_pop( offspring )
        , m_alleles( alleles )
        , m_parents( first, last )
        , m_offspring_index( off_idx )
        , m_recomb_rate( recomb_rate )
        , m_seq_bias( seq_bias )
    {}

    batch_crossover_task( const self_type & other ) :
        m_rng( other.m_rng )
        , m_parent_pop( other.m_parent_pop )
        , m_offspring_pop( other.m_offspring_pop )
        , m_alleles( other.m_alleles )
        , m_parents( other.m_parents )
        , m_offspring_index( other.m_offspring_index )
        , m_recomb_rate( other.m_recomb_rate )
        , m_seq_bias( other.m_seq_bias )
    { }

    void operator()() {

        event_distribution_type     event_dist( m_recomb_rate);
        boost::random::bernoulli_distribution< double > bias_dist( m_seq_bias);

        const_iterator mate_it = m_parents.begin(), mate_end = m_parents.end();

        const unsigned int PARENT_ALLELE_COUNT = m_parent_pop->getMaxAlleles();
        const unsigned int PARENT_SIZE = m_parent_pop->haploid_genome_count();

        const unsigned int OFFSPRING_ALLELE_COUNT = m_offspring_pop->getMaxAlleles();
        const unsigned int OFFSPRING_SIZE = m_offspring_pop->haploid_genome_count();

        size_t i = m_offspring_index;

        while( mate_it != mate_end ) {

            individual_pointer p = m_parent_pop->getIndividual( mate_it->first );
            genome_pointer c = m_offspring_pop->getHaploidGenome( i++ );

            event_type evts;
            fill_events( evts, event_dist( *m_rng ) );
            classifier_type cfier0( &m_alleles->getPositions(), evts );
            bool _swap = bias_dist( *m_rng );

            run_crossover_task( cfier0, p, PARENT_ALLELE_COUNT, PARENT_SIZE, c, OFFSPRING_ALLELE_COUNT, OFFSPRING_SIZE, _swap );

            p = m_parent_pop->getIndividual( mate_it->second );
            c = m_offspring_pop->getHaploidGenome( i++ );

            event_type evts1;
            fill_events( evts1, event_dist( *m_rng ) );
            classifier_type cfier1( &m_alleles->getPositions(), evts1 );
            _swap = bias_dist( *m_rng );

            run_crossover_task( cfier1, p, PARENT_ALLELE_COUNT, PARENT_SIZE, c, OFFSPRING_ALLELE_COUNT, OFFSPRING_SIZE, _swap );

            ++mate_it;
        }
    }

    virtual ~batch_crossover_task() {}

protected:

    inline void fill_events( event_type & evt, unsigned int N ) {
        while( N-- ) {
            evt.push_back( m_pos_dist( *m_rng ) );
        }
    }

    event_type make_events( unsigned int N ) {
        event_type res;

        while( N-- ) {
            res.push_back( m_pos_dist( *m_rng ) );
        }

        return res;
    }

    void run_crossover_task( const classifier_type & cls, individual_pointer parent, const unsigned int PARENT_ALLELE_COUNT, const unsigned int PARENT_STEP, genome_pointer offspring, const unsigned int OFFSPRING_ALLELE_COUNT, const unsigned int OFFSPRING_STEP, bool should_swap_strands ) {
        if( cls.event_count() == 0 ) {
            if( should_swap_strands ) {
                ++parent;
            }

            unsigned int i = 0;

            while( i < PARENT_ALLELE_COUNT)  {
                *offspring = *parent;
                parent += PARENT_STEP;
                offspring += OFFSPRING_STEP;
                i += bit_helper_type::BITS_PER_BLOCK;
            }

            // fill 
//            while( i < OFFSPRING_ALLELE_COUNT) {
//                *offspring = bit_helper_type::ALL_UNSET;
//                offspring += OFFSPRING_STEP;
//                i += bit_helper_type::BITS_PER_BLOCK;
//            }
        } else if( should_swap_strands ) {
            crossover_task_swapped( cls, parent, PARENT_ALLELE_COUNT, PARENT_STEP, offspring, OFFSPRING_ALLELE_COUNT, OFFSPRING_STEP);
        } else {
            crossover_task( cls, parent, PARENT_ALLELE_COUNT, PARENT_STEP, offspring, OFFSPRING_ALLELE_COUNT, OFFSPRING_STEP);
        }
    }

    void crossover_task( const classifier_type & cls, individual_pointer parent, const unsigned int PARENT_ALLELE_COUNT, const unsigned int PARENT_STEP, genome_pointer offspring, const unsigned int OFFSPRING_ALLELE_COUNT, const unsigned int OFFSPRING_STEP ) {

//        std::cerr << "Parent Population Bounds: " << PARENT_ALLELE_COUNT << " alleles x " << PARENT_STEP << " genomes" << std::endl;
//        std::cerr << "Offspring Population Bounds: " << OFFSPRING_ALLELE_COUNT << " alleles x " << OFFSPRING_STEP << " genomes" << std::endl;
        crossover_type xover( cls );

        unsigned int i = 0;
        while( i < PARENT_ALLELE_COUNT ) {
            const block_type t = *parent;
            const block_type b = *(parent + 1);

            //std::cerr << "checking allele range: [" << i << ", " << i + bit_helper_type::BITS_PER_BLOCK << ")" << std::endl; 
            *offspring = xover.crossover( t, b, i );               

            parent += PARENT_STEP;
            offspring += OFFSPRING_STEP;
            i += bit_helper_type::BITS_PER_BLOCK;
        }

//        while( i < OFFSPRING_ALLELE_COUNT ) {
//            *offspring = bit_helper_type::ALL_UNSET;
//            offspring += OFFSPRING_STEP;
//            i += bit_helper_type::BITS_PER_BLOCK;
//        }
    }

    void crossover_task_swapped( const classifier_type & cls, individual_pointer parent, const unsigned int PARENT_ALLELE_COUNT, const unsigned int PARENT_STEP, genome_pointer offspring, const unsigned int OFFSPRING_ALLELE_COUNT, const unsigned int OFFSPRING_STEP ) {
        crossover_type xover( cls );

        unsigned int i = 0;
        while( i < PARENT_ALLELE_COUNT) {
            const block_type b = *parent;
            const block_type t = *(parent + 1);

            *offspring = xover.crossover( t, b, i );               

            parent += PARENT_STEP;
            offspring += OFFSPRING_STEP;
            i += bit_helper_type::BITS_PER_BLOCK;
        }

//        while( i < OFFSPRING_ALLELE_COUNT ) {
//            *offspring = bit_helper_type::BITS_PER_BLOCK;
//            offspring += OFFSPRING_STEP;
//            i += bit_helper_type::BITS_PER_BLOCK;
//        }
    }

    random_engine_type  * m_rng;
    space_type * m_parent_pop, * m_offspring_pop;
    allele_type         * m_alleles;

    mate_pair_type m_parents;

    unsigned int m_offspring_index;

    double m_recomb_rate, m_seq_bias;

    position_distribution_type  m_pos_dist;
};

}   // namespace genetics
}   // namesapce clotho

#endif  // CLOTHO_BATCH_CROSSOVER_TASK_POPULATION_SPACE_COLUMNAR_HPP_
