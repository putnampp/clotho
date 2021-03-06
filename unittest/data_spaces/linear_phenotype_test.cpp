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
#include <boost/test/unit_test.hpp>

#include "../unittest_config.h"

#include "clotho/data_spaces/allele_space/allele_space.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"
#include "clotho/data_spaces/phenotype_evaluator/linear_combination_method.hpp"
#include "clotho/data_spaces/phenotype_evaluator/phenotype_evaluator.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/property_tree/ptree.hpp>
#include <vector>
#include <unordered_map>

typedef boost::random::mt19937                                                  random_engine_type;
typedef double                                                                  position_type;
typedef double                                                                  weight_type;
typedef unsigned long long                                                      block_type;
typedef clotho::genetics::qtl_allele_vectorized< position_type, weight_type >   allele_type;
typedef clotho::genetics::genetic_space< allele_type, block_type, clotho::genetics::column_aligned >              genetic_space_type;
typedef clotho::genetics::TraitWeightAccumulator< genetic_space_type >          trait_accum_type;

typedef typename genetic_space_type::individual_genome_type                     individual_type;

typedef double                                                                      phenotype_type;
typedef clotho::genetics::linear_combination< trait_accum_type, phenotype_type >    linear_combination_type;
typedef clotho::genetics::phenotype_evaluator< linear_combination_type >            phenotype_eval_type;

typedef clotho::genetics::AlleleGenerator< random_engine_type, allele_type >   allele_generator_type;

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( single_trait_linear_combination_phenotype_test ) {
    
    size_t exp_alleles = 200, exp_genomes = 20, exp_traits = 1;
    genetic_space_type gs;

    gs.grow( exp_genomes, exp_alleles );
    gs.getAlleleSpace().grow( exp_alleles, exp_traits );
    gs.getSequenceSpace().clear();

    random_engine_type          rand(1234);
    boost::property_tree::ptree exp_config;

    exp_config.put( "neutral.p", 0.5 );     // qtl alleles may be neutral
    exp_config.put( "traits.mean", 0.0 );
    exp_config.put( "traits.sigma", 1.0 );

    allele_generator_type   agen( &rand, exp_config );

    size_t i = 0;
    while( i < exp_alleles ) {
        agen( gs.getAlleleSpace(), i++ );
    }

    size_t id = 10;
    individual_type exp_ind = std::make_pair(2 * id, 2 * id + 1);
    individual_type obs_ind = gs.getIndividualAt( id );

    BOOST_REQUIRE_MESSAGE( obs_ind.first == exp_ind.first, "Unexpected individual's first sequence; Observed: " << obs_ind.first << "; Expected: " << exp_ind.first );
    BOOST_REQUIRE_MESSAGE( obs_ind.second == exp_ind.second, "Unexpected individual's second sequence; Observed: " << obs_ind.second << "; Expected: " << exp_ind.second );

    typedef  std::vector< std::pair< size_t, size_t> > allele_index_type;
    typedef typename allele_index_type::iterator allele_iterator;

    std::unordered_map< size_t, weight_type > exp_trait_w;
    allele_index_type exp_a;

    exp_a.push_back( std::make_pair( exp_ind.first, 12) );
    exp_a.push_back( std::make_pair( exp_ind.first, 35) );
    exp_a.push_back( std::make_pair( exp_ind.first, 120) );
    exp_a.push_back( std::make_pair( exp_ind.second, 50) );
    exp_a.push_back( std::make_pair( exp_ind.second, 150) );
    exp_a.push_back( std::make_pair( exp_ind.second, 199) );

    phenotype_type exp_pheno = 0.0;

    allele_iterator ait = exp_a.begin();
    while( ait != exp_a.end() ) {
        std::pair< size_t, size_t > tmp = *ait++;
        gs.getSequenceSpace().flip( tmp.first, tmp.second );

        typedef typename allele_type::trait_iterator iterator;

        iterator it = gs.getAlleleSpace().getTraitIterator( tmp.second );
        weight_type w = 0.0;
        while( it.hasNext() ) {
            weight_type _w = it.next();
            w += _w;
        }

        exp_pheno += w;

        if( exp_trait_w.find( tmp.first ) != exp_trait_w.end() ) {
            w += exp_trait_w[ tmp.first ];
        }
        exp_trait_w[ tmp.first ] = w;
    }

    trait_accum_type tacc;

    tacc.update( gs );

    phenotype_eval_type ph;

    ph.update( &gs, tacc );

    weight_type exp_weight_0 = exp_trait_w[ exp_ind.first ];
    weight_type exp_weight_1 = exp_trait_w[ exp_ind.second ];

    weight_type obs_weight_0 = tacc.getTraitAt( exp_ind.first )->at( 0 );
    weight_type obs_weight_1 = tacc.getTraitAt( exp_ind.second )->at( 0 );

    BOOST_REQUIRE_MESSAGE( exp_weight_0 == obs_weight_0, "Unexpected trait weight 0; Observed: " << obs_weight_0 << "; Expected: " << exp_weight_0 );
    BOOST_REQUIRE_MESSAGE( exp_weight_1 == obs_weight_1, "Unexpected trait weight 0; Observed: " << obs_weight_1 << "; Expected: " << exp_weight_1 );

    phenotype_type obs_pheno = ph.getPhenotypeAt( id );

    BOOST_REQUIRE_MESSAGE( exp_pheno != 0.0, "Unexpected zero phenotype value" );
    BOOST_REQUIRE_MESSAGE( exp_pheno == obs_pheno, "Unexpected phenotype; Observed: " << obs_pheno << "; Expected: " << exp_pheno );
}

BOOST_AUTO_TEST_SUITE_END()
