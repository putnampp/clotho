#ifndef INFINITE_SITE_GENERATOR_HPP_
#define INFINITE_SITE_GENERATOR_HPP_

template < class Sequence >
struct infinite_site {
    typedef Sequence type;
};

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/powerset/vector_subset.hpp"
#include "clotho/powerset/powerset_no_dup_pred.hpp"
#include "clotho/utility/random_generator.hpp"

namespace clotho {
namespace utility {

template < class URNG, class E, class B, class BM, class EK >
class random_generator< URNG, infinite_site< clotho::powersets::variable_subset< E, B, BM, EK > > > {
public:
    typedef E result_type;
    typedef clotho::powersets::variable_subset< E, B, BM, EK > sequence_type;

    typedef clotho::mutations::no_duplicate_pred< typename sequence_type::powerset_type > no_dup_type;
    typedef clotho::utility::random_generator< URNG, E >    random_element;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_gen( rng, config ) {
    }

    result_type operator()( std::shared_ptr< sequence_type > & seq, unsigned int age = 0 ) {
        return operator()( *seq, age );
    }

    result_type operator()( const sequence_type & seq, unsigned int age = 0 ) {
        result_type res = m_gen(age);

        no_dup_type tester( seq.getParent() );
        while( !tester( res ) ) {
            res = m_gen();
        }

        return res;
    }

protected:
    URNG * m_rng;
    random_element m_gen;
};

template < class URNG, class E, class B, class BM, class EK >
class random_generator< URNG, infinite_site< clotho::powersets::vector_subset< E, B, BM, EK > > > {
public:
    typedef E result_type;
    typedef clotho::powersets::vector_subset< E, B, BM, EK > sequence_type;

    typedef clotho::mutations::no_duplicate_pred< typename sequence_type::powerset_type > no_dup_type;
    typedef clotho::utility::random_generator< URNG, E >    random_element;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_gen( rng, config ) {
    }

    result_type operator()( std::shared_ptr< sequence_type > & seq, unsigned int age = 0 ) {
        return operator()( *seq, age );
    }

    result_type operator()( const sequence_type & seq, unsigned int age = 0 ) {
        result_type res = m_gen( age );

        no_dup_type tester( seq.getParent() );
        while( !tester( res ) ) {
            res = m_gen();
        }

        return res;
    }

protected:
    URNG * m_rng;
    random_element m_gen;
};
}   // namespace utility
}   // namespace clotho

#endif  // INFINITE_SITE_HPP_
