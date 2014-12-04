#ifndef SEQUENCE_HELPER_HPP_
#define SEQUENCE_HELPER_HPP_

#include "clotho/utility/iterator_helper.hpp"

template < class IndividualType >
struct sequence_helper {
    typedef typename IndividualType::sequence_type type;
    typedef clotho::utility::iterator_helper< IndividualType > iterator_helper_type;
};

template < class SequenceType >
struct sequence_helper< std::pair< SequenceType, SequenceType > > {
    typedef SequenceType        sequence_type;
    typedef clotho::utility::iterator_helper< std::pair< SequenceType, SequenceType > > iterator_helper_type;
};

#endif  // SEQUENCE_HELPER_HPP_
