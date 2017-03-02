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
#ifndef OFFSPRING_GENERATOR_POPULATION_SPACE_INDEX_VECTOR_ALIGNMENT_HPP_
#define OFFSPRING_GENERATOR_POPULATION_SPACE_INDEX_VECTOR_ALIGNMENT_HPP_

#include "clotho/data_spaces/offspring_generator/offspring_generator_def.hpp"
#include "clotho/data_spaces/population_space/population_spaces.hpp"

#include "clotho/data_spaces/crossover/common_crossovers.hpp"
#include "clotho/data_spaces/phenotype_evaluator/common_phenotype_accumulators.hpp"

#include "clotho/data_spaces/selection/selection_generator_factory.hpp"

#include "clotho/utility/bit_helper.hpp"

#include <map>

namespace clotho {
namespace genetics {

template < class RNG, class IndexType, class BlockType, class AlleleSpaceType, class FitnessType, class TraitSpaceType, class FreeBufferType >
class offspring_generator<RNG, population_space< index_vector_alignment< IndexType, BlockType >, TraitSpaceType >, AlleleSpaceType, FitnessType, TraitSpaceType, FreeBufferType > {
public:
    typedef offspring_generator< RNG, population_space< index_vector_alignment< IndexType, BlockType >, TraitSpaceType >, AlleleSpaceType, FitnessType, TraitSpaceType, FreeBufferType > self_type;

    typedef RNG                                             random_engine_type;
    typedef index_vector_alignment< IndexType, BlockType >  alignment_type;
    typedef TraitSpaceType                                  trait_space_type;

    typedef population_space< alignment_type, trait_space_type >   population_type;

    typedef typename alignment_type::index_type             index_type;
    typedef typename alignment_type::block_type             block_type;

    typedef AlleleSpaceType                                 allele_type;

    typedef std::vector< unsigned int >                     mutation_pool_type;
    typedef std::vector< unsigned int >                     mutation_distribution_type;

    typedef typename population_type::genome_pointer        genome_pointer;
    typedef typename population_type::row_pointer           row_pointer;

    typedef FreeBufferType                                  free_buffer_type;

    typedef FitnessType                                     fitness_type;

    typedef clotho::genetics::common_crossover< RNG, population_type, allele_type > crossover_type;
    typedef clotho::genetics::common_phenotype_accumulator< population_type, trait_space_type > phenotype_method;

    typedef typename phenotype_method::weight_type          weight_type;

    typedef clotho::utility::BitHelper< block_type >        bit_helper_type;

    typedef std::vector< std::pair< unsigned long long, unsigned long long > >               time_vector;

    offspring_generator( random_engine_type * rng, allele_type * alleles, mutation_pool_type * mut_pool, mutation_distribution_type * mut_dist, fitness_type * fit, trait_space_type * traits, double recomb_rate, double bias_rate ) :
        m_rng( rng )
        , m_parents( NULL )
        , m_offspring( NULL )
        , m_alleles( alleles )
        , m_mut_pool( mut_pool )
        , m_mut_dist( mut_dist )
        , m_fitness( fit )
        , m_off_begin( 0 )
        , m_off_end( 0 )
        , m_all_neutral( 0 )
        , m_crossover_method( rng, alleles, recomb_rate, bias_rate )
        , m_pheno_method( traits )
    {}

    void reset( population_type * parents, population_type * offspring, const free_buffer_type & buf, unsigned int off_idx, unsigned int off_end, bool allNeutral ) {
        m_parents = parents;
        m_offspring = offspring;

        m_free_space.reset( buf );

        m_off_begin = off_idx;
        m_off_end = off_end;

        m_all_neutral = allNeutral;
    }

    void operator()() {
        // crossover parents of offspring population
        //
        timer_type xover_time;
        crossover();
        xover_time.stop();
        
        // mutate offspring population
        //
        timer_type mutate_time;
        mutate();
        mutate_time.stop();

        //
        // evaluate phenotype of  offspring population
        //
        timer_type pheno_time;
        phenotype();
        pheno_time.stop();

        // evaluate fixed alleles with offspring population
        //
        timer_type fixed_time;
        fixed();
        fixed_time.stop();

        // persist time logs
        m_crossover_times.push_back( std::make_pair( xover_time.getStart(), xover_time.getStop()));
        m_mutate_times.push_back( std::make_pair( mutate_time.getStart(), mutate_time.getStop()));
        m_phenotype_times.push_back( std::make_pair( pheno_time.getStart(), pheno_time.getStop()));
        m_fixed_times.push_back( std::make_pair( fixed_time.getStart(), fixed_time.getStop()));
    }

    void record( boost::property_tree::ptree & xo, boost::property_tree::ptree & mt, boost::property_tree::ptree & ph, boost::property_tree::ptree & fx ) {

        boost::property_tree::ptree xo_begins, xo_ends;
        buildTimeLog( xo_begins, xo_ends, m_crossover_times );
        xo.put_child( "start", xo_begins);
        xo.put_child( "stop", xo_ends);

        boost::property_tree::ptree mt_begins, mt_ends;
        buildTimeLog( mt_begins, mt_ends, m_mutate_times );
        mt.put_child( "start", mt_begins);
        mt.put_child( "stop", mt_ends);

        boost::property_tree::ptree ph_begins, ph_ends;
        buildTimeLog( ph_begins, ph_ends, m_phenotype_times );
        ph.put_child( "start", ph_begins);
        ph.put_child( "stop", ph_ends);

        boost::property_tree::ptree fx_begins, fx_ends;
        buildTimeLog( fx_begins, fx_ends, m_fixed_times );
        fx.put_child( "start", fx_begins);
        fx.put_child( "stop", fx_ends);
    }

    virtual ~offspring_generator() {}

protected:

    void buildTimeLog( boost::property_tree::ptree & begins, boost::property_tree::ptree & ends, time_vector & t ) {
        for( typename time_vector::iterator it = t.begin(); it != t.end(); it++ ) {
            clotho::utility::add_value_array( begins, it->first );
            clotho::utility::add_value_array( ends, it->second );
        }
    }

    void crossover() {
        std::shared_ptr< selection_details< random_engine_type, unsigned int > > sel = clotho::genetics::make_selection_generator( m_rng, m_fitness);

        for( unsigned int i = m_off_begin; i < m_off_end; ++i ) {
            typename selection_details< random_engine_type >::parent_pair res = (*sel)();

            perform_crossover_method( 2 * res.first, 2 * i );
            perform_crossover_method( 2 * res.second, 2 * i + 1 );
        }
    }

    void perform_crossover_method( unsigned int p_idx, unsigned int o_idx ) {
        genome_pointer top = m_parents->begin_genome( p_idx ), top_end = m_parents->end_genome(p_idx);
        genome_pointer bottom = m_parents->begin_genome( p_idx + 1 ), bottom_end = m_parents->end_genome( p_idx + 1 );
        row_pointer off = m_offspring->begin_row( o_idx );

        off->clear();

        m_crossover_method( top, top_end, bottom, bottom_end, off ); 
    }

    void mutate() {
        for( unsigned int i = 2 * m_off_begin; i < 2 * m_off_end; ++i ) {
            unsigned int lo = m_mut_dist->at( i ), hi = m_mut_dist->at( i + 1 );

            while( lo < hi ) {
                unsigned int all_idx = m_mut_pool->at( lo );
                m_offspring->mutate( i, all_idx );
                ++lo;
            }
        } 
    }

    void phenotype() {
        if( !m_all_neutral ) {
            for( unsigned int i = 2 * m_off_begin; i < 2 * m_off_end; ++i ) {
                genome_pointer first = m_offspring->begin_genome( i ), last = m_offspring->end_genome( i );

                m_pheno_method( first, last, m_alleles );

                m_offspring->updateGenomeWeights(i, m_pheno_method.getResults());
            }
        } else {
            // constant phenotype
            m_pheno_method.resetBuffer();

            for( unsigned int i = 2 * m_off_begin; i < 2 * m_off_end; ++i )
                m_offspring->updateGenomeWeights( i, m_pheno_method.getResults() ); 
        }
    }

    void fixed() {
        std::map< typename population_type::index_type, unsigned int > idx_set;

        for( unsigned int i = 2 * m_off_begin; i < 2 * m_off_end; ++i ) {
            genome_pointer first = m_offspring->begin_genome(i), last = m_offspring->end_genome(i);

            while( first != last ) {

                typename std::map< typename population_type::index_type, unsigned int >::iterator it = idx_set.find(*first);
                if( it == idx_set.end() ) {
                    idx_set.insert( std::make_pair( *first, 1 ) );
                } else {
                    it->second++;
                }

                ++first;
            }

        }

        unsigned int M = m_offspring->getMaxBlocks();
        block_type * tmp = new block_type[ 2 * M ];
        block_type * fix_tmp = tmp + M;
        memset( tmp, 0, 2 * sizeof( block_type ) * M );
        for( typename std::map< typename population_type::index_type, unsigned int >::iterator it = idx_set.begin(); it != idx_set.end(); ++it ) {
            unsigned int idx = ( it->first / bit_helper_type::BITS_PER_BLOCK );
            block_type m = ((block_type) 1 << ( it->first % bit_helper_type::BITS_PER_BLOCK ) );

            tmp[ idx ] |= m;
            if( it->second == 2 * (m_off_end - m_off_begin)) {
                fix_tmp[ idx ] |= m;
            }
        }

        m_free_space.update( tmp, M );
        m_free_space.update( fix_tmp, M );

        delete [] tmp;
    }

    random_engine_type * m_rng;

    population_type * m_parents, * m_offspring;

    allele_type * m_alleles;

    mutation_pool_type * m_mut_pool;
    mutation_distribution_type * m_mut_dist;

    fitness_type  * m_fitness;

    unsigned int m_off_begin, m_off_end;

    bool m_all_neutral;

    crossover_type m_crossover_method;

    phenotype_method m_pheno_method;
    
    free_buffer_type  m_free_space;

    time_vector m_crossover_times, m_mutate_times, m_phenotype_times, m_fitness_times, m_fixed_times;
};

}   // namespace genetics
}   // namespace clotho

#endif  // OFFSPRING_GENERATOR_POPULATION_SPACE_INDEX_VECTOR_ALIGNMENT_HPP_
