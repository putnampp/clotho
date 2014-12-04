#ifndef SEQUENCE_TRAIT_ACCUMULATOR_HPP_
#define SEQUENCE_TRAIT_ACCUMULATOR_HPP_

#include "trait_accumulator.hpp"
#include "trait_helper.hpp"
#include "clotho/utility/iterator_helper.hpp"

template < class SequenceType >
class sequence_trait_accumulator {
public:
    typedef sequence_trait_accumulator< SequenceType >          self_type;
    typedef SequenceType                                        sequence_type;

    typedef typename trait_helper< sequence_type >::element_type element_type;
    typedef typename trait_helper< sequence_type >::weight_type weight_type;

    typedef trait_accumulator< element_type, weight_type >      trait_accumulator_type;
    typedef typename trait_accumulator_type::result_type        result_type;

    typedef clotho::utility::iterator_helper< sequence_type >   elem_helper;
    typedef typename elem_helper::type                          elem_iterator;

    sequence_trait_accumulator( result_type & r) : t_acc(r) {}
    sequence_trait_accumulator( const self_type & other ) : t_acc( other.t_acc ) {}

    void operator()( sequence_type & seq ) {
        elem_iterator first = elem_helper::make_first(seq), last = elem_helper::make_last(seq);

        std::for_each( first, last, t_acc );
    }

protected:
    trait_accumulator_type t_acc;
};

template < class SequenceType >
class sequence_trait_accumulator< std::shared_ptr< SequenceType > > : 
    public sequence_trait_accumulator< SequenceType > {
public:
    typedef sequence_trait_accumulator< std::shared_ptr< SequenceType > > self_type;
    typedef sequence_trait_accumulator< SequenceType >  base_accumulator;
    typedef std::shared_ptr< SequenceType >             sequence_type;

    typedef typename base_accumulator::result_type      result_type;

    sequence_trait_accumulator( result_type & r ) : base_accumulator( r ) {}
    sequence_trait_accumulator( const self_type & other ) : base_accumulator( other ) {}

    void operator()( sequence_type & seq ) {
        base_accumulator::operator()( *seq );
    }
};

#endif  // SEQUENCE_TRAIT_ACCUMULATOR_HPP_
