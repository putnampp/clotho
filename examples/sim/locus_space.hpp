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
#ifndef LOCUS_SPACE_HPP_
#define LOCUS_SPACE_HPP_

#include "space.hpp"
#include "locus.hpp"

class locus_space : public space {
public:

    locus_space();

    virtual size_t size();
    virtual void prune();

    virtual ~locus_space();

protected:
    std::vector< std::shared_ptr< locus > > m_loci;
};

#endif  // LOCUS_SPACE_HPP_
