#ifndef DISCRETE_SELECTOR_HPP_
#define DISCRETE_SELECTOR_HPP_

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

namespace clotho {
namespace selection {

template < class IndividualPointer, class EnvironmentPointer >
class DiscreteSelection {
public:
    DiscreteSelection( gsl_rng * r, double * fitnesses, size_t s ) :
        m_rng( r ),
        m_lookup( NULL ) {
        m_lookup = gsl_ran_discrete_preproc( s, fitnesses );
    }

    size_t operator()() {
        return gsl_ran_discrete( m_rng, m_lookup );
    }

    //std::pair< individual_pointer , individual_pointer > operator()( environment_type * env, double f = 0.0 ) {
    std::pair< IndividualPointer, IndividualPointer > operator()( EnvironmentPointer env, double f = 0.0 ) {
        size_t i0 = gsl_ran_discrete( m_rng, m_lookup );
        size_t i1 = ((gsl_rng_uniform(m_rng) <= f ) ? i0 : gsl_ran_discrete( m_rng, m_lookup ));

        std::pair< IndividualPointer, IndividualPointer > res = std::make_pair( env->at(i0), env->at(i1));

        return res;
    }

    virtual ~DiscreteSelection() {
        gsl_ran_discrete_free( m_lookup );
    }
protected:
    gsl_rng * m_rng;
    gsl_ran_discrete_t * m_lookup;
};

}   // namespace selection
}   // namespace clotho

#endif  // DISCRETE_SELECTOR_HPP_
