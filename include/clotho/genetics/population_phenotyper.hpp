#ifndef POPULATION_PHENOTYPER_HPP_
#define POPULATION_PHENOTYPER_HPP_

template < class Population >
class population_phenotyper;

#include <vector>
#include <unordered_map>

#include "clotho/genetics/individual_trait_accumulator.hpp"
#include "clotho/genetics/sequence_helper.hpp"

template < class IndividualType >
class population_phenotyper< std::vector< IndividualType > > {
public:

    typedef population_phenotyper< std::vector< IndividualType > > self_type;

    typedef individual_trait_accumulator< IndividualType >      phenotype_accum_type;
    typedef typename phenotype_accum_type::result_type          phenotype_type;

    typedef typename phenotype_accum_type::sequence_accumulator_type    sequence_effect;
    typedef typename sequence_effect::result_type                       effect_type;
    
    typedef typename phenotype_accum_type::sequence_helper_type sequence_helper_type;
    typedef typename sequence_helper_type::sequence_type        sequence_type;
    typedef typename sequence_helper_type::iterator_helper_type seq_iterator_helper;
    typedef typename seq_iterator_helper::range_type        seq_iterator_range;

    typedef std::unordered_map< sequence_type, effect_type > sequence_effect_map;
    typedef typename sequence_effect_map::iterator           effect_iterator;

    typedef std::vector< phenotype_type >           result_type;

    population_phenotyper( ) : m_pop_res( NULL ) {}
    population_phenotyper( result_type & pop_res ) : m_pop_res(&pop_res) {}

    population_phenotyper( const self_type & other ) : m_pop_res( other.m_pop_res ) {}

    phenotype_type operator()( IndividualType & ind ) {
    //    seq_iterator first = seq_iterator_helper::make_first(ind), last = seq_iterator_helper::make_last(ind);

        seq_iterator_range ran = seq_iterator_helper::make_range( ind );
        phenotype_type                  pheno_res;
        //phenotype_accum_type            ind_accum(pheno_res);

        while( ran.first != ran.second ) {
            effect_iterator it = m_effects.find( *ran.first );
            if( it == m_effects.end() ) {
                effect_type tmp_effect;
                sequence_effect seff( tmp_effect );

                seff( *ran.first );

                std::pair< effect_iterator, bool > tmp = m_effects.insert( std::make_pair( *ran.first, tmp_effect ));

                it = tmp.first;
            }

            //ind_accum( it->second );
            typename phenotype_type::iterator pit = pheno_res.begin();
            typename phenotype_type::iterator first = it->second.begin();
            while( first != it->second.end() ) {
                if( pit == pheno_res.end() ) {
                    while( first != it->second.end() ) {
                        pheno_res.push_back( *first );
                        ++first;
                    }
                    break;
                }
                (*pit) += *first;
                ++pit;
                ++first;
            }
            ++ran.first;
        }

        if( m_pop_res ) m_pop_res->push_back( pheno_res );

        return pheno_res;
    }

protected:
    result_type         * m_pop_res;
    sequence_effect_map m_effects;
};

#endif  // POPULATION_PHENOTYPER_HPP_
