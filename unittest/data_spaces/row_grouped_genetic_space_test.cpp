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

#include "clotho/data_spaces/allele_space/qtl_allele.hpp"
#include "clotho/data_spaces/population_space/genetic_space.hpp"
#include "clotho/utility/bit_helper.hpp"

BOOST_AUTO_TEST_SUITE( test_data_space )

BOOST_AUTO_TEST_CASE( row_aligned_genetic_space_size_test ) {
    typedef clotho::genetics::qtl_allele_vectorized< double, double >   allele_type;
    typedef unsigned long long                                          block_type;
    typedef clotho::genetics::genetic_space< allele_type, block_type, clotho::genetics::row_grouped< 2 > >  genetic_space_type;

    size_t _alleles = 10, _genomes = 11;
    size_t _seqs = 2 * _genomes;
    genetic_space_type gs;

    gs.grow( _genomes, _alleles );

    size_t obs_alleles = gs.allele_count();
    size_t obs_genomes = gs.genome_count();
    size_t obs_seqs = gs.sequence_count();

    BOOST_REQUIRE_MESSAGE( obs_alleles == _alleles, "Unexpected number of alleles in genetic space; Observed: " << obs_alleles << "; Expected: " << _alleles );
    BOOST_REQUIRE_MESSAGE( obs_genomes == _genomes, "Unexpected number of genomes in genetic space; Observed: " << obs_genomes << "; Expected: " << _genomes );
    BOOST_REQUIRE_MESSAGE( obs_seqs == _seqs, "Unexpected number of sequences in genetic space; Observed: " << obs_seqs << "; Expected: " << _seqs );

    size_t exp_alloc_size = _seqs * (_alleles / clotho::utility::BitHelper< block_type >::BITS_PER_BLOCK + 1);
    size_t obs_alloc_size = gs.getSequenceSpace().allocated_size();
    
    BOOST_REQUIRE_MESSAGE( exp_alloc_size == obs_alloc_size, "Unexpected allocated size for genetic sequence space; Observed: " << obs_alloc_size << "; Expected: " << exp_alloc_size );

    size_t exp_traits = 1;
    size_t obs_traits = gs.getAlleleSpace().trait_count();

    BOOST_REQUIRE_MESSAGE( exp_traits == obs_traits, "Unexpected number of traits allocated; Observed: " << obs_traits << "; Expected: " << exp_traits );

}

BOOST_AUTO_TEST_SUITE_END()
