#ifndef INDIVIDUAL_REPRODUCTION_HPP_
#define INDIVIDUAL_REPRODUCTION_HPP_

template < class IndividualType, class MutationGenerator, class RecombinationGenerator >
class individual_reproduction {
public:
    typedef IndividualType individual_type;
    typedef IndividualType result_type;

    result_type operator()( individual_type & p0, individual_type & p1 ) {
        assert(false);
    }
};

#include <utility>

template < class SequenceType, class MutationGenerator, class RecombinationGenerator >
class individual_reproduction< std::pair< SequenceType, SequenceType >, MutationGenerator, RecombinationGenerator > {
public:
    typedef SequenceType    sequence_type;
    typedef std::pair< SequenceType, SequenceType > individual_type;
    typedef individual_type result_type;

    typedef MutationGenerator                           mutation_generator_type;
    typedef RecombinationGenerator                      recombination_generator_type;

    typedef typename MutationGenerator::result_type        mutation_type;
    typedef typename RecombinationGenerator::result_type   recombination_type;

    individual_reproduction( mutation_generator_type & mut, recombination_generator_type & rec ) :
        m_mut_gen( mut )
        , m_rec_gen( rec ) {
    }

    result_type operator()( individual_type & p0, individual_type & p1, unsigned int gen = 0 ) {
        mutation_type   mut0 = m_mut_gen( gen );
        recombination_type rec0 = m_rec_gen();

        sequence_type   c0 = rec0( p0, mut0.event_count() == 0 );
        mut0( c0 );

        mutation_type mut1 = m_mut_gen( gen );
        recombination_type rec1 = m_rec_gen();

        sequence_type   c1 = rec1( p1, mut1.event_count() == 0 );
        mut1(c1);

        return std::make_pair(c0, c1);
    }
protected:
    mutation_generator_type         m_mut_gen;
    recombination_generator_type    m_rec_gen;
};

#endif  // INDIVIDUAL_REPRODUCTION_HPP_
