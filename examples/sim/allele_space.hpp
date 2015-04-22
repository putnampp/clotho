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
#ifndef ALLELE_SPACE_HPP_
#define ALLELE_SPACE_HPP_

#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "allele.hpp"

/**
 * Allele Space maps implicitly to Effect Space.  That is, the allele
 * at index i of AlleleSpace maps to the effect_size at index i in
 * EffectSpace.
 *
 * Alleles maps explicitly to loci in LocusSpace.
 *
 */
class AlleleSpace {
public:
    AlleleSpace( boost::property_tree::ptree & config );

    std::shared_ptr< Allele > create();

    size_t size();
    void prune();

    virtual ~AlleleSpace();
protected:
    EffectSpace                                 m_effects;
    LocusSpace                                  m_loci;
    std::vector< std::shared_ptr< Allele > >    m_alleles;

};

#endif  // ALLELE_SPACE_HPP_
