#ifndef INDIVIDUAL_GENERATOR_HPP_
#define INDIVIDUAL_GENERATOR_HPP_

template < class Population, class SelectionModel, class ReproductionModel >
class individual_generator;

#include "individual_selector.hpp"
#include "individual_reproduction.hpp"

template < class IndividualType, class URNG, class MutationModel, class RecombinationModel >
class individual_generator< std::vector< IndividualType >, individual_selector< URNG >, individual_reproduction< IndividualType, MutationModel, RecombinationModel > > {
public:
    typedef std::vector< IndividualType >   population_type;
    typedef IndividualType                  individual_type;
    typedef IndividualType                  result_type;

    typedef individual_selector< URNG >     selection_type;
    typedef individual_reproduction< IndividualType, MutationModel, RecombinationModel > reproduction_type;

    individual_generator( population_type * pop, selection_type & sel, reproduction_type & repro ) :
        m_pop( pop )
        , m_sel( sel )
        , m_repro( repro ) {
    }

    result_type operator()() {
        individual_type p0 = m_pop->at(m_sel()), p1 = m_pop->at( m_sel());

        return m_repro(p0, p1);
    }
protected:
    population_type * m_pop;
    selection_type  m_sel;
    reproduction_type m_repro;
};


#endif  // INDIVIDUAL_GENERATOR_HPP_
