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
#ifndef INDIVIDUAL_HPP_
#define INDIVIDUAL_HPP_

#include <vector>
#include <memory>

#include "sequence.hpp"

struct individual_param {
    unsigned char nParents;
    unsigned char nPloidy;
};

/**
 * individual interface
 *
 * All individuals are expected to have at least one parent
 * All individuals are expected to have at least one genetic sequence
 *
 * Not expecting individuals to have children, basically because if you
 * have an entire population and you know their parents, then you can figure
 * out the children/siblings.
 */
struct individual {

    /**
     * pointer to the parent from which sequence(s) are inheritted
     */
    virtual std::shared_ptr< individual > get_parent( unsigned int ) = 0;

    /**
     * pointer to the sequence which was inheritted from a parent.
     */
    virtual std::shared_ptr< sequence > get_sequence( unsigned int ) = 0;
};

class Individual : public individual {
public:

    Individual( const individual_param & params );

    virtual std::shared_ptr< individual >   get_parent( unsigned int idx );
    virtual std::shared_ptr< sequence >     get_sequence( unsigned int idx );

    virtual ~Individual();
protected:
    std::vector< std::shared_ptr< individual > > m_parents;
    std::vector< std::shared_ptr< sequence > > m_genetics;
};

#endif  // INDIVIDUAL_HPP_
