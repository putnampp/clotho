#ifndef INDIVIDUAL_REPRODUCTION_DEF_HPP_
#define INDIVIDUAL_REPRODUCTION_DEF_HPP_

template < class IndividualType, class MutationGenerator, class RecombinationGenerator, class Tag >
class individual_reproduction {
public:
    typedef IndividualType individual_type;
    typedef IndividualType result_type;

    result_type operator()( individual_type & p0, individual_type & p1 ) {
        assert(false);
    }
};

#endif  // INDIVIDUAL_REPRODUCTION_DEF_HPP_
