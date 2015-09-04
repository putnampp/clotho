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
#ifndef POPULATION_SPACE_HPP_
#define POPULATION_SPACE_HPP_

template < class SSpaceType, class ASpaceType >
struct PopulationSpace {
    typedef PopulationSpace< SSpaceType, ASpaceType > self_type;
    typedef SSpaceType  sequence_space_type;
    typedef ASpaceType  allele_space_type;

    sequence_space_type sequences;
    allele_space_type   alleles;

    template < class EventSpaceType >
    void resize( self_type * pop, EventSpaceType * mut_events, unsigned int seqs ) {
        alleles.expand_relative_to( pop->alleles, mut_events );

        sequences.resize( alleles.get_device_space(), seqs );
    }
};

#endif  // POPULATION_SPACE_HPP_
