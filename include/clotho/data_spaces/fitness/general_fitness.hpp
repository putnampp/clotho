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
#ifndef CLOTHO_GENERAL_FITNESS_HPP_
#define CLOTHO_GENERAL_FITNESS_HPP_

#include <boost/property_tree/ptree.hpp>

#include "clotho/fitness/fitness_toolkit.hpp"

#include "clotho/data_spaces/population_space/population_spaces.hpp"
#include "clotho/data_spaces/phenotype_evaluator/batch_phenotype_mt.hpp"

#include <vector>

namespace clotho {
namespace genetics {

class GeneralFitness {
public:
    typedef std::shared_ptr< ifitness_generator >           fitness_generator;
    typedef std::shared_ptr< ifitness >                     fitness_operator;

    typedef unsigned int                                    individual_id_type;

    typedef typename ifitness::result_type                  fitness_type;
    typedef std::vector< fitness_type >                     fitness_vector;

    typedef typename fitness_vector::iterator               iterator;
    typedef typename fitness_vector::const_iterator         const_iterator;

    GeneralFitness( boost::property_tree::ptree & config ) :
        m_fit_gen()
    {

        m_fit_gen = fitness_toolkit::getInstance()->get_tool( config );

        if( !m_fit_gen ) {
            fitness_toolkit::getInstance()->tool_configurations( config );
        }
    }

    template < class PhenotypeSpaceType >
    void operator()( const PhenotypeSpaceType & phenos ) {
        typedef typename PhenotypeSpaceType::phenotype_type         weight_type;

        const unsigned int N = phenos.individual_count();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        weight_type * tmp = phenos.getPhenotypes();
        unsigned int traits = phenos.trait_count();

        unsigned int i = 0;
        while( i < N ) {
            fitness_type score = (*op)(tmp, tmp + traits);

            setFitness( i, score );

            tmp += traits;
            ++i;
        }
    }

    template < class BlockType, class WeightType, class TraitSpaceType >
    void operator()( BatchPhenotypeMT< population_space_columnar< BlockType, WeightType >, TraitSpaceType > & phenos ) {
        typedef BatchPhenotypeMT< population_space_columnar< BlockType, WeightType >, TraitSpaceType > phenotype_space_type;
        const unsigned int N = phenos.individual_count();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        for( unsigned int i = 0; i < N; ++i ) {
            typename phenotype_space_type::phenotype_iterator first = phenos.begin_individual_phenotype(i), last = phenos.end_individual_phenotype(i);

            typename phenotype_space_type::phenotype_vector w( first, last );

            fitness_type score = (*op)(w);

            setFitness( i, score );
        }
    }

    template < class BlockType, class WeightType, class TraitSpaceType >
    void operator()( BatchPhenotypeMT< population_space_row< BlockType, WeightType >, TraitSpaceType > & phenos ) {
        typedef BatchPhenotypeMT< population_space_row< BlockType, WeightType >, TraitSpaceType > phenotype_space_type;
        const unsigned int N = phenos.individual_count();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        for( unsigned int i = 0; i < N; ++i ) {
            typename phenotype_space_type::phenotype_iterator first = phenos.begin_individual_phenotype(i), last = phenos.end_individual_phenotype(i);

            typename phenotype_space_type::phenotype_vector w( first, last );

            fitness_type score = (*op)(w);

            setFitness( i, score );
        }
    }

    template < class BlockType, class WeightType >
    void operator()( population_space_row< BlockType, WeightType > * pop ) {
        const unsigned int N = pop->getIndividualCount();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        typedef typename population_space_row< BlockType, WeightType >::const_weight_iterator iterator;
        for( unsigned int i = 0; i < N; ++i ) {
            iterator f = pop->begin_genome_traits( 2 * i ), e = pop->end_genome_traits( 2 * i );
            typename population_space_row< BlockType, WeightType >::weight_vector w( f, e );

            f = pop->begin_genome_traits( 2 * i + 1);
            e = pop->end_genome_traits( 2 * i + 1 );

            unsigned int j = 0;
            while( f != e ) {
                w[ j++ ] += *f++;
            }

            fitness_type score = (*op)(w);

            setFitness( i, score );
        }
    }

    template < class BlockType, class WeightType >
    void operator()( population_space< row_block_alignment< BlockType >, trait_space_vector< WeightType > > * pop ) {
        const unsigned int N = pop->getIndividualCount();

        if( N == 0 ) return;

        fitness_operator op = m_fit_gen->generate( N );

        typedef typename population_space< row_block_alignment< BlockType >, trait_space_vector< WeightType > >::const_weight_iterator iterator;
        for( unsigned int i = 0; i < N; ++i ) {
            iterator f = pop->begin_genome_traits( 2 * i ), e = pop->end_genome_traits( 2 * i );
            typename population_space_row< BlockType, WeightType >::weight_vector w( f, e );

            f = pop->begin_genome_traits( 2 * i + 1);
            e = pop->end_genome_traits( 2 * i + 1 );

            unsigned int j = 0;
            while( f != e ) {
                w[ j++ ] += *f++;
            }

            fitness_type score = (*op)(w);

            setFitness( i, score );
        }
    }

    template < class Iterator >
    void operator()( Iterator first, Iterator last ) {
        const unsigned int N = last - first;

        fitness_operator op = m_fit_gen->generate( N );

        unsigned int i = 0;
        while( first != last ) {
            fitness_type score = (*op)( *first );

            setFitness( i++, score );
            ++first;
        }
    }

    unsigned int individual_count() const {
        return m_fitness.size();
    }

    void setFitness( unsigned int i, fitness_type score ) {
        if( i < m_fitness.size() ) {
            m_fitness[ i ] = score;
        } else {
            do {
                m_fitness.push_back( score );
            } while( m_fitness.size() <= i );
        }
    }

    fitness_type getFitness( unsigned int idx ) {
#ifdef DEBUGGING
        assert( idx < m_fitness.size() ) ;
#endif  // DEBUGGING
        return m_fitness[ idx ];
    }

    iterator begin() {
        return m_fitness.begin();
    }

    iterator end() {
        return m_fitness.end();
    }

    const_iterator begin() const {
        return m_fitness.begin();
    }

    const_iterator end() const {
        return m_fitness.end();
    }

    virtual ~GeneralFitness() {}

protected:
    fitness_generator   m_fit_gen;

    fitness_vector      m_fitness;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_GENERAL_FITNESS_HPP_
