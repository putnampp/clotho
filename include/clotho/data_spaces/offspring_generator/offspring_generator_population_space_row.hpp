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
#ifndef OFFSPRING_GENERATOR_POPULATION_SPACE_ROW_HPP_
#define OFFSPRING_GENERATOR_POPULATION_SPACE_ROW_HPP_

#include "clotho/data_spaces/offspring_generator/offspring_generator_def.hpp"
#include "clotho/data_spaces/population_space/population_space_row.hpp"

#include "clotho/data_spaces/crossover/common_crossover.hpp"
#include "clotho/data_spaces/phenotype_evaluator/common_phenotype_accumulator.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class BlockType, class WeightType, class AlleleSpaceType, class SelectionType, class TraitSpaceType, class FreeBufferType >
class offspring_generator<RNG, population_space_row< BlockType, WeightType >, AlleleSpaceType, SelectionType, TraitSpaceType, FreeBufferType > {
public:
    typedef offspring_generator< RNG, population_space_row< BlockType, WeightType >, AlleleSpaceType, SelectionType, TraitSpaceType, FreeBufferType > self_type;

    typedef RNG                                             random_engine_type;
    typedef population_space_row< BlockType, WeightType >   population_type;
    typedef AlleleSpaceType                                 allele_type;

    typedef std::vector< unsigned int >                     mutation_pool_type;
    typedef std::vector< unsigned int >                     mutation_distribution_type;

    typedef typename population_type::block_type            block_type;
    typedef typename population_type::genome_pointer        genome_pointer;

    typedef FreeBufferType                                  free_buffer_type;

    typedef SelectionType                                   selection_type;
    typedef TraitSpaceType                                  trait_space_type;

    typedef clotho::genetics::common_crossover< RNG, population_type, allele_type > crossover_type;

    typedef clotho::genetics::common_phenotype_accumulator< population_type, trait_space_type > phenotype_method;

    typedef typename phenotype_method::weight_type          weight_type;

    typedef std::vector< std::pair< unsigned long long, unsigned long long > >               time_vector;

    offspring_generator( random_engine_type * rng, population_type * parents, population_type * offspring, allele_type * alleles, mutation_pool_type * mut_pool, mutation_distribution_type * mut_dist, selection_type * sel, trait_space_type * traits, const free_buffer_type & fbuf, unsigned int off_idx, unsigned int off_end, double recomb_rate, double bias_rate, bool allNeutral ) :
        m_rng( rng )
        , m_parents( parents )
        , m_offspring( offspring )
        , m_alleles( alleles )
        , m_mut_pool( mut_pool )
        , m_mut_dist( mut_dist )
        , m_select( sel )
        , m_off_begin( off_idx )
        , m_off_end( off_end )
        , m_all_neutral( allNeutral )
        , m_crossover_method( rng, alleles, recomb_rate, bias_rate )
        , m_pheno_method( traits, alleles->getNeutrals() )
        , m_free_space( fbuf )
    {
    }

    void operator()() {
        // crossover parents of offspring population
        //
        m_crossover_timer.start();
        crossover();
        m_crossover_timer.stop();
        
        // mutate offspring population
        //
        m_mutate_timer.start();
        mutate();
        m_mutate_timer.stop();
        //
        // evaluate phenotype of  offspring population
        //
        m_phenotype_timer.start();
        phenotype();
        m_phenotype_timer.stop();

        // evaluate fixed alleles with offspring population
        //
        m_fixed_timer.start();
        fixed();
        m_fixed_timer.stop();
    }

    void record( time_vector & xover, time_vector & mut, time_vector & pheno, time_vector & fixed ) {
        xover.push_back( std::make_pair( m_crossover_timer.getStart(), m_crossover_timer.getStop() ) );
        mut.push_back( std::make_pair( m_mutate_timer.getStart(), m_mutate_timer.getStop() ) );
        pheno.push_back( std::make_pair( m_phenotype_timer.getStart(), m_phenotype_timer.getStop() ) );
        fixed.push_back( std::make_pair( m_fixed_timer.getStart(), m_fixed_timer.getStop() ) );
    }

    void recordCrossover( boost::property_tree::ptree & log ) {
        recordTime( log, m_crossover_timer );
    }

    void recordMutate( boost::property_tree::ptree & log ) {
        recordTime( log, m_mutate_timer );
    }

    void recordPhenotype( boost::property_tree::ptree & log ) {
        recordTime( log, m_phenotype_timer );
    }

    void recordFixed( boost::property_tree::ptree & log ) {
        recordTime( log, m_fixed_timer );
    }

    virtual ~offspring_generator() {}

protected:
    void recordTime( boost::property_tree::ptree & log, const timer_type & t ) {
        boost::property_tree::ptree elapsed, start, stop;

        elapsed = log.get_child("elapsed", elapsed);
        start = log.get_child( "start", start );
        stop = log.get_child( "stop", stop );

        clotho::utility::add_value_array( elapsed, t );
        clotho::utility::add_value_array( start, t.getStart() );
        clotho::utility::add_value_array( stop, t.getStop() );
        
        log.put_child( "elapsed", elapsed );
        log.put_child( "start", start );
        log.put_child( "stop", stop );
    }

    void crossover() {
        typename selection_type::iterator b = m_select->begin() + m_off_begin, e = m_select->begin() + m_off_end;

        unsigned int off_idx = 2 * m_off_begin;
        while( b != e ) {
            perform_crossover_method( 2 * b->first, off_idx++ );
            perform_crossover_method( 2 * b->second, off_idx++ );

            ++b;
        }
    }

    void perform_crossover_method( unsigned int p_idx, unsigned int o_idx ) {
        genome_pointer top = m_parents->begin_genome( p_idx ), top_end = m_parents->end_genome(p_idx);
        genome_pointer bottom = m_parents->begin_genome( p_idx + 1 ), bottom_end = m_parents->end_genome( p_idx + 1 );
        genome_pointer off = m_offspring->begin_genome( o_idx ), off_end = m_offspring->end_genome(o_idx);

        m_crossover_method( top, top_end, bottom, bottom_end, off, off_end ); 
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

                m_pheno_method( first, last );

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
        unsigned int M = m_offspring->getMaxBlocks();

        for( unsigned int i = 2 * m_off_begin; i < 2 * m_off_end; ++i ) {
            block_type * f = m_offspring->begin_genome( i );
            m_free_space.update( f, M );
        }
    }

    random_engine_type * m_rng;

    population_type * m_parents, * m_offspring;

    allele_type * m_alleles;

    mutation_pool_type * m_mut_pool;
    mutation_distribution_type * m_mut_dist;

    selection_type  * m_select;

    const unsigned int m_off_begin, m_off_end;

    bool m_all_neutral;

    crossover_type m_crossover_method;

    phenotype_method m_pheno_method;
    
    free_buffer_type m_free_space;

    timer_type  m_crossover_timer, m_mutate_timer, m_phenotype_timer, m_fitness_timer, m_fixed_timer;
};

}   // namespace genetics
}   // namespace clotho

#endif  // OFFSPRING_GENERATOR_POPULATION_SPACE_ROW_HPP_
