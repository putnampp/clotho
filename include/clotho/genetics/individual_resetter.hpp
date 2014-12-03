#ifndef INDIVIDUAL_RESETTER_HPP_
#define INDIVIDUAL_RESETTER_HPP_

template < class IndividualType >
struct individual_resetter;

#include "resetter.hpp"

template < class SequenceType >
struct individual_resetter< std::pair< SequenceType, SequenceType > > {
    typedef std::pair< SequenceType, SequenceType > individual_type;
    typedef resetter< SequenceType > resetter_type;

    void operator()( individual_type & ind ) {
        resetter_type::reset( ind.first );
        resetter_type::reset( ind.second  );;
    }
};

#endif  // INDIVIDUAL_RESETTER_HPP_
