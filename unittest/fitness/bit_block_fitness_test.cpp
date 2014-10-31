#include <boost/test/unit_test.hpp>

#include "../unittest_config.h"

#include "clotho/fitness/bit_block_fitness.hpp"

#include <vector>
#include <iterator>

typedef unsigned long   Block;
typedef double          Result;

template < class EIterator >
void buildElements( EIterator eit ) {
    for( unsigned int i = 0; i < sizeof( Block ) * 8; ++i) {
        (*eit++) = (Result) i;
    }
}

struct basic_het {
    void operator()( Result & res, Result v ) {
        res += v;
    }
};

BOOST_AUTO_TEST_SUITE( test_fitness )

BOOST_AUTO_TEST_CASE( no_fit_test ) {
    typedef clotho::fitness::no_fit fit_algo_type;

    typedef clotho::fitness::bit_block_fitness< fit_algo_type, fit_algo_type, fit_algo_type, Result > fitness_type;

    fitness_type fit;

    std::vector< double > elems;
    buildElements( std::back_inserter(elems));

    BOOST_REQUIRE_MESSAGE( elems.size() == sizeof(Block) * 8, "Unexpected number of elements");

    double f = 1.0;
    double res = fit( f, (Block) 0xA, (Block) 0xB, (Block) 8, elems.begin() );

    BOOST_REQUIRE_MESSAGE( res == f, "Unexpected result returned");
}

BOOST_AUTO_TEST_CASE( just_het ) {
    typedef clotho::fitness::no_fit fit_algo_type;

    typedef clotho::fitness::bit_block_fitness< basic_het, fit_algo_type, fit_algo_type, Result > fitness_type;

    fitness_type fit;

    std::vector< double > elems;
    buildElements( std::back_inserter(elems) );
    double f = 0.;
    double f_exp = (0. + 0.);   // 0xA ^ 0xB = 0x1
    double res = fit( f, (Block) 0xA, (Block) 0xB, (Block) 8, elems.begin() );

    BOOST_REQUIRE_MESSAGE( res == f_exp, "Unexpected result returned");
}

BOOST_AUTO_TEST_SUITE_END()
