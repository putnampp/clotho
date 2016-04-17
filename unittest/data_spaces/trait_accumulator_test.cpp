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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include "clotho/data_spaces/allele_space/qtl_allele.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/data_spaces/phenotype_evaluator/trait_accumulator.hpp"

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( trait_accum_size_test ) {

    typedef clotho::genetics::qtl_allele_vectorized< double, double >           allele_type;
    typedef unsigned long long                                                  block_type;
    typedef clotho::genetics::genetic_space< allele_type, block_type >          genetic_space_type;
    typedef clotho::genetics::TraitWeightAccumulator< genetic_space_type >      trait_accumulator_type;

    size_t exp_alleles = 12, exp_genomes = 110, exp_traits = 15;
    genetic_space_type  gs;

    gs.grow( exp_genomes, exp_alleles );
    gs.getAlleleSpace().grow( exp_alleles, exp_traits );
    
    trait_accumulator_type weights( gs );

    size_t obs_size = weights.size();

    BOOST_REQUIRE_MESSAGE( obs_size == exp_genomes * 2, "Unexpected number of weights; Observed: " << obs_size << "; Expected: " << (2 * exp_genomes));

    size_t i = 0;

    while( i < obs_size ) {
        size_t obs_traits = weights.getTraitAt(i).size();

        BOOST_REQUIRE_MESSAGE( obs_traits == exp_traits, "Unexpected number of traits at index " << i << "; Observed: " << obs_traits << "; Expected: " << exp_traits );

        ++i;
    }
}

BOOST_AUTO_TEST_CASE( trait_accum_update_test ) {

    typedef double                                                              weight_type;
    typedef clotho::genetics::qtl_allele_vectorized< double, weight_type >      allele_type;
    typedef unsigned long long                                                  block_type;
    typedef clotho::genetics::genetic_space< allele_type, block_type >          genetic_space_type;
    typedef clotho::genetics::TraitWeightAccumulator< genetic_space_type >      trait_accumulator_type;

    typedef typename allele_type::trait_iterator                                trait_iterator;

    typedef typename trait_accumulator_type::trait_vector_type                  trait_vector_type;
    typedef typename trait_accumulator_type::trait_helper_type                  trait_helper_type;

    size_t exp_alleles = 12, exp_genomes = 110, exp_traits = 15;
    genetic_space_type  gs;

    gs.grow( exp_genomes, exp_alleles );
    gs.getAlleleSpace().grow( exp_alleles, exp_traits );

    gs.getSequenceSpace().clear();

    boost::random::mt19937  rand_engine;
    boost::random::uniform_01< weight_type >  uni;

    size_t i = 0, j = 0;


    // flip some bits
    while( i < 2 * exp_genomes ) {

        gs.getSequenceSpace().flip( i, (j % exp_alleles));

        j += (i % exp_alleles);
        ++i;
    }

    i = 0;
    while( i < exp_alleles ) {

        j = 0;
        while( j < exp_traits ) {
            weight_type w = uni( rand_engine );
            gs.getAlleleSpace().updateTraitWeight( i, j, w );
            ++j;
        }
        ++i;
    } 

    trait_accumulator_type weights( gs );

    size_t obs_size = weights.size();

    BOOST_REQUIRE_MESSAGE( obs_size == exp_genomes * 2, "Unexpected number of weights; Observed: " << obs_size << "; Expected: " << (2 * exp_genomes));

    i = 0;
    while( i < obs_size ) {
        size_t obs_traits = weights.getTraitAt(i).size();

        BOOST_REQUIRE_MESSAGE( obs_traits == exp_traits, "Unexpected number of traits at index " << i << "; Observed: " << obs_traits << "; Expected: " << exp_traits );

        ++i;
    }

    i = 0;
    while ( i < exp_genomes * 2 ) {

        j = 0;
        trait_vector_type res = trait_helper_type::makeEmptyTrait( exp_traits );
        while( j < exp_alleles ) {
            
            if( gs.getSequenceSpace()( i, j ) ) {
                trait_iterator it = gs.getAlleleSpace().getTraitIterator( j );
                size_t k = 0;
                while( it.hasNext() ) {
                    weight_type r = it.next();
                    res[ k++ ] += r;
                }
            }
            ++j;
        }

        bool eq = true;
        j = 0;
        while( eq && (j < exp_traits )) {
            eq = (res[ j ] == weights.getTraitAt(i)[j]);
            ++j;
        }

        BOOST_REQUIRE_MESSAGE( eq, "Unexpected weight vector at " << i);
        
        i += 11;
    }
}

BOOST_AUTO_TEST_SUITE_END()
