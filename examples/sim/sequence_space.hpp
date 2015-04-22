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
#ifndef SEQUENCE_SPACE_HPP_
#define SEQUENCE_SPACE_HPP_

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/dynamic_bitset.hpp>

#include "allele_space.hpp"

class SequenceSpace {
public:
    typedef boost::dynamic_bitset< unsigned long >  sequence_type;
    typedef std::shared_ptr< sequence_type >        sequence_ptr;

    SequenceSpace( boost::property_tree::ptree & config );

    sequence_ptr create();

    sequence_ptr clone( sequence_ptr s );

    size_t size();
    void prune();

    virtual ~SequenceSpace();
protected:
    AlleleSpace m_alleles;

    std::vector< sequence_ptr > m_seq_map;
};

#endif  // SEQUENCE_SPACE_HPP_
