#ifndef INDIVIDUAL_GENERATOR_HPP_
#define INDIVIDUAL_GENERATOR_HPP_

template < class Population, class SelectionModel, class ReproductionModel >
class individual_generator;

#include "individual_selector.hpp"
#include "individual_reproduction.hpp"

template < class IndividualType, class URNG, class MutationModel, class RecombinationModel, class ReproMethodTag >
class individual_generator< std::vector< IndividualType >, individual_selector< URNG >, individual_reproduction< IndividualType, MutationModel, RecombinationModel, ReproMethodTag > > {
public:
    typedef std::vector< IndividualType >   population_type;
    typedef IndividualType                  individual_type;
    typedef IndividualType                  result_type;

    typedef individual_selector< URNG >     selection_type;
    typedef individual_reproduction< IndividualType, MutationModel, RecombinationModel, ReproMethodTag > reproduction_type;

    individual_generator( population_type * pop, selection_type & sel, reproduction_type & repro, unsigned int age = 0 ) :
        m_pop( pop )
        , m_sel( sel )
        , m_repro( repro )
        , m_gen( age ) {
    }

    result_type operator()() {
        individual_type p0 = m_pop->at(m_sel()), p1 = m_pop->at( m_sel());

        return m_repro(p0, p1, m_gen);
    }
protected:
    population_type * m_pop;
    selection_type  m_sel;
    reproduction_type m_repro;
    unsigned int    m_gen;
};


#endif  // INDIVIDUAL_GENERATOR_HPP_
