//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef MUTATE_RECOMBINE_INDIVIDUAL_REPRODUCTION_HPP_
#define MUTATE_RECOMBINE_INDIVIDUAL_REPRODUCTION_HPP_

#include "clotho/genetics/individual_reproduction_def.hpp"

struct mutate_recombine_tag {};

template < class SequenceType, class MutationGenerator, class RecombinationGenerator >
class individual_reproduction< std::pair< SequenceType, SequenceType >, MutationGenerator, RecombinationGenerator, mutate_recombine_tag > {
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
        // TODO: Very 'inefficient' method; used for proof of concept
        // Inefficient because it generates (and adds mutations) which may not be passed along to child as a result of recombination
        //
        // assumes sequence_type is shared_ptr; hence copy amounts to returning sequence
        // introduce new mutations into either of parent 0's sequences
        mutation_type   mut = m_mut_gen( gen );
        sequence_type mp00 = p0.first;
        if( mut.event_count() > 0 ) {
            mp00 = p0.first->clone();
            mut(mp00);
        }

        mut = m_mut_gen(gen);
        sequence_type mp01 = p0.second;
        if(mut.event_count() > 0) {
            mp01 = p0.second->clone();
            mut(mp01);
        }

        recombination_type rec = m_rec_gen();

        // since sequences have already been cloned when mutated, should always copy sequence if recombination results in an input sequence
        sequence_type   c0 = rec( std::make_pair( mp00, mp01 ), true );

        // introduce new mutations into either of parent 1's sequences
        mut = m_mut_gen(gen);
        sequence_type mp10 = p1.first;
        if(mut.event_count() > 0) {
            mp10 = p1.first->clone();
            mut(mp10);
        }

        mut = m_mut_gen(gen);
        sequence_type mp11 = p1.second;
        if( mut.event_count() > 0 ) {
            mp11 = p1.second->clone();
            mut( mp11 );
        }

        rec = m_rec_gen();
        sequence_type   c1 = rec( std::make_pair( mp10, mp11), true );

        return std::make_pair(c0, c1);
    }

    virtual ~individual_reproduction() {}
protected:
    mutation_generator_type         m_mut_gen;
    recombination_generator_type    m_rec_gen;
};
#endif  // MUTATE_RECOMBINE_INDIVIDUAL_REPRODUCTION_HPP_
