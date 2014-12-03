#ifndef SIMULATE_ENGINE_BASE_HPP_
#define SIMULATE_ENGINE_BASE_HPP_

#include "clotho/powerset/variable_subset.hpp"
#include "clotho/utility/random_generator.hpp"

#include "clotho/classifiers/region_classifier.hpp"

#include <vector>
#include <utility>

template < class URNG, class AlleleType, class LogType, class TimerType >
class simulate_engine_base {
public:
    typedef URNG                                                        rng_type;
    typedef AlleleType                                                  allele_type;
    typedef LogType                                                     log_type;
    typedef TimerType                                                   timer_type;

    typedef clotho::utility::random_generator< rng_type, allele_type >  allele_generator;

    typedef clotho::powersets::variable_subset< allele_type, unsigned long > sequence_type;
    typedef typename sequence_type::pointer                                  sequence_pointer;
    typedef typename sequence_type::powerset_type                            allele_set_type;

    typedef clotho::classifiers::region_classifier< allele_type >            classifier_type;
    typedef clotho::utility::random_generator< rng_type, classifier_type >   classifier_generator;

    typedef std::pair< sequence_pointer, sequence_pointer >                  individual_type;
    typedef std::vector< individual_type >                                   population_type;
    typedef typename population_type::iterator                               population_iterator;
};

#endif  // SIMULATE_ENGINE_BASE_HPP_
