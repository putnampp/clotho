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
#ifndef ALLELE_HPP_
#define ALLELE_HPP_

#include <memory>

#include "locus.hpp"

class Allele {
public:
    typedef unsigned int    allele_state_type;

    Allele();

    std::shared_ptr< locus >        getLocus();
    std::shared_ptr< effect_size >  getEffectSize();

    virtual ~Allele();

protected:
    std::shared_ptr< locus >        m_locus;
    std::shared_ptr< effect_size >  m_effect;
    allele_state_type               m_state;
};

#endif  // ALLELE_HPP_
