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
#ifndef CLOTHO_POPULATION_INITIALIZER_HPP_
#define CLOTHO_POPULATION_INITIALIZER_HPP_

#include <boost/property_tree/ptree.hpp>

#include <vector>

namespace clotho {
namespace genetics {

class PopulationInitializer {
public:

    PopulationInitializer( boost::property_tree::ptree & config );

    bool hasDistribution() const;
    size_t getAlleleCount() const;

    std::vector< double >::iterator begin();
    std::vector< double >::const_iterator begin() const;

    std::vector< double >::iterator end();
    std::vector< double >::const_iterator end() const;

    virtual ~PopulationInitializer();

private:
    std::vector< double > m_allele_freq;
};

} // namespace genetics
} // namespace clotho

#endif  // CLOTHO_POPULATION_INITIALIZER_HPP_
