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
#ifndef BATCH_OFFSPRING_GENERATOR_HPP_
#define BATCH_OFFSPRING_GENERATOR_HPP_

#include "clotho/data_spaces/crossover/common_crossovers.hpp"
#include "clotho/data_spaces/phenotype_evaluator/common_phenotype_accumulators.hpp"

#include "clotho/data_spaces/selection/selection_generator_factory.hpp"

#include "clotho/recombination/sequence_bias_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

namespace clotho {
namespace genetics {

template < class RNG, class PopulationType, class AlleleSpaceType, class FitnessType, class TraitSpaceType, class FreeBufferType >
class batch_offspring_generator {
public:

    typedef RNG                     random_engine_type;
    typedef PopulationType          population_type;

    typedef AlleleSpaceType         allele_type;

    typedef std::vector< unsigned int >                     mutation_pool_type;
    typedef std::vector< unsigned int >                     mutation_distribution_type;

    typedef FitnessType             fitness_type;

    struct batch_parameters {
        population_type * parent, * offspring;
        fitness_type * parent_fit, * off_fit;

        mutation_pool_type * mut_pool;
        mutation_distribution_type * mut_dist;

        unsigned int p_start, p_end;
        unsigned int o_start, o_end;

        bool all_neutral;

        batch_parameters( population_type * p
            , population_type * o
            , fitness_type * p_fit
            , fitness_type * o_fit
            , mutation_pool_type * pool
            , mutation_distribution_type * dist
            , unsigned int p_s
            , unsigned int p_e
            , unsigned int o_s
            , unsigned int o_e
            , bool neutral ) :
            parent( p )
            , offspring( o )
            , parent_fit( p_fit)
            , off_fit( o_fit )
            , mut_pool( pool )
            , mut_dist( dist )
            , p_start( p_s )
            , p_end( p_e )
            , o_start(o_s)
            , o_end(o_e)
            , all_neutral( neutral )
        {}

        batch_parameters( const batch_parameters & bp ) :
            parent(bp.parent)
            , offspring( bp.offspring )
            , parent_fit( bp.parent_fit )
            , off_fit( bp.off_fit )
            , mut_pool( bp.mut_pool )
            , mut_dist( bp.mut_dist )
            , p_start( bp.p_start )
            , p_end( bp.p_end )
            , o_start( bp.o_start )
            , o_end( bp.o_end )
            , all_neutral( bp.all_neutral )
        {}

        virtual ~batch_parameters() {}
    };

    typedef FreeBufferType          free_buffer_type;
    typedef TraitSpaceType          trait_space_type;

    typedef clotho::genetics::common_crossover< RNG, population_type, allele_type > crossover_type;
    typedef clotho::genetics::common_phenotype_accumulator< population_type, trait_space_type > phenotype_method;
    typedef typename phenotype_method::weight_type          weight_type;

    typedef batch_parameters                                parameter_type;
    typedef std::vector< parameter_type * >                 parameter_vector;

    typedef std::vector< std::pair< unsigned long long, unsigned long long > >               time_vector;

    typedef typename population_type::genome_pointer        genome_pointer;

    batch_offspring_generator( random_engine_type * rng, allele_type * alleles, trait_space_type * traits, boost::property_tree::ptree & config ) :
        m_rng( rng )
        , m_alleles( alleles )
        , m_recomb_rate( config )
        , m_bias_rate( config )
        , m_crossover_method( rng, alleles, m_recomb_rate.m_rho, m_bias_rate.m_bias )
        , m_pheno_method( traits )
    {}

    void reset( const free_buffer_type & buf ) {
        m_pop_params.clear();
        m_free_space.reset( buf );
    }

    void addPopulation( parameter_type * p ) {
        m_pop_params.push_back( p );
    }

    void operator()() {
        if( m_pop_params.empty() ) return;

        for( unsigned int i = 0; i < m_pop_params.size(); ++i ) {

            parameter_type * params = m_pop_params[ i ];

            std::shared_ptr< selection_details< random_engine_type, unsigned int> > sel = buildSelection(params);

            timer_type xover_time;
            crossover( params, sel );
            xover_time.stop();

            timer_type mutate_time;
            mutate( params );
            mutate_time.stop();

            timer_type pheno_time;
            phenotype( params );
            pheno_time.stop();

            timer_type fit_time;
            fitness( params );
            fit_time.stop();

            m_crossover_times.push_back( std::make_pair( xover_time.getStart(), xover_time.getStop()));
            m_mutate_times.push_back( std::make_pair( mutate_time.getStart(), mutate_time.getStop()));
            m_phenotype_times.push_back( std::make_pair( pheno_time.getStart(), pheno_time.getStop()));
            m_fitness_times.push_back( std::make_pair( fit_time.getStart(), fit_time.getStop() ) );
            m_fixed_times.push_back( std::make_pair( fit_time.getStop() + 5, fit_time.getStop() + 10));
        }
        
        parameter_type * param = m_pop_params[ m_pop_params.size() - 1 ];
        
        timer_type fixed_time;
        fixed( param );
        fixed_time.stop();

        m_fixed_times.push_back( std::make_pair( fixed_time.getStart(), fixed_time.getStop()));
    }

    void record( boost::property_tree::ptree & xo, boost::property_tree::ptree & mt, boost::property_tree::ptree & ph, boost::property_tree::ptree & fx, boost::property_tree::ptree & ft ) {

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

        boost::property_tree::ptree ft_begins, ft_ends;
        buildTimeLog( ft_begins, ft_ends, m_fitness_times );
        ft.put_child( "start", ft_begins);
        ft.put_child( "stop", ft_ends);
    }

    virtual ~batch_offspring_generator() {}

protected:

    std::shared_ptr< selection_details< random_engine_type, unsigned int> > buildSelection( parameter_type * params ) {
        assert( params->parent_fit != NULL );

        std::shared_ptr< selection_details< random_engine_type, unsigned int> > sel = clotho::genetics::make_selection_generator( m_rng, params->parent_fit, params->p_start, params->p_end );
        return sel;
    }

    inline void buildTimeLog( boost::property_tree::ptree & begins, boost::property_tree::ptree & ends, time_vector & t ) {
        for( typename time_vector::iterator it = t.begin(); it != t.end(); it++ ) {
            clotho::utility::add_value_array( begins, it->first );
            clotho::utility::add_value_array( ends, it->second );
        }
    }

    void crossover( parameter_type * params, std::shared_ptr< selection_details< random_engine_type, unsigned int> > sel ) {
        for( unsigned int i = params->o_start; i < params->o_end; ++i ) {
            typename selection_details< random_engine_type >::parent_pair res = (*sel)();

            perform_crossover_method( params->parent, params->offspring, 2 * (res.first + params->p_start), 2 * i );
            perform_crossover_method( params->parent, params->offspring, 2 * (res.second + params->p_start), 2 * i + 1 );
        }

    }

    void perform_crossover_method( population_type * parent, population_type * offspring, unsigned int p_idx, unsigned int o_idx ) {
        genome_pointer top = parent->begin_genome( p_idx ), top_end = parent->end_genome( p_idx );
        genome_pointer bottom = parent->begin_genome( p_idx + 1), bottom_end = parent->end_genome( p_idx + 1 );
        genome_pointer off = offspring->begin_genome(o_idx), off_end = offspring->end_genome( o_idx );

        m_crossover_method( top, top_end, bottom, bottom_end, off, off_end );
    }

    void mutate( parameter_type * params ) {
        for( unsigned int i = 2 * params->o_start; i < 2 * params->o_end; ++i ) {
            unsigned int lo = params->mut_dist->at( i ), hi = params->mut_dist->at( i + 1 );

            while( lo < hi ) {
                unsigned int all_idx = params->mut_pool->at( lo );
                params->offspring->mutate( i, all_idx );
                ++lo;
            }
        }
    }

    void phenotype( parameter_type * params ) {
        if( !params->all_neutral ) {
            for( unsigned int i = 2 * params->o_start; i < 2 * params->o_end; ++i ) {
                genome_pointer first = params->offspring->begin_genome( i ), last = params->offspring->end_genome( i );

                m_pheno_method( first, last, m_alleles->getNeutrals() );

                params->offspring->updateGenomeWeights(i, m_pheno_method.getResults());
            }
        } else {
            // constant phenotype
            m_pheno_method.resetBuffer();

            for( unsigned int i = 2 * params->o_start; i < 2 * params->o_end; ++i )
                params->offspring->updateGenomeWeights( i, m_pheno_method.getResults() );
        }
    }

    void fitness( parameter_type * params ) {
        (*params->off_fit)(params->offspring, params->o_start, params->o_end );
    }

    void fixed( parameter_type * params ) {

        unsigned int M = params->offspring->getMaxBlocks();
        unsigned int lo = 2 * params->o_start, hi = 2 * params->o_end;
        while( lo < hi ) {
            genome_pointer f = params->offspring->begin_genome( lo++ );
            m_free_space.update( f, M );
        }
    }

    random_engine_type  * m_rng;

    allele_type         * m_alleles;

    parameter_vector    m_pop_params;

    recombination_rate_parameter< double > m_recomb_rate;
    sequence_bias_parameter< double > m_bias_rate;

    crossover_type      m_crossover_method;
    phenotype_method    m_pheno_method;

    free_buffer_type    m_free_space;

    time_vector m_crossover_times, m_mutate_times, m_phenotype_times, m_fitness_times, m_fixed_times;
};

}   // namespace genetics
}   // namespace clotho

#endif  // BATCH_OFFSPRING_GENERATOR_HPP_

