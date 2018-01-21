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

#include "clotho/initializer/PopulationInitializer.hpp"

namespace clotho {
namespace genetics {

PopulationInitializer::PopulationInitializer( boost::property_tree::ptree & config ) {
    boost::property_tree::ptree  init_tree, all_dist;
    init_tree = config.get_child( "initialize", init_tree);

    all_dist = init_tree.get_child( "allele_distribution", all_dist );
    for (auto& item : all_dist ) {
        m_allele_freq.push_back( item.second.get_value<double>() );
    }
}

bool PopulationInitializer::hasDistribution() const {
    return !m_allele_freq.empty();
}

size_t PopulationInitializer::getAlleleCount() const {
    return m_allele_freq.size();
}

std::vector< double >::iterator PopulationInitializer::begin() {
    return m_allele_freq.begin();
}

std::vector< double >::const_iterator PopulationInitializer::begin() const {
    return m_allele_freq.begin();
}

std::vector< double >::iterator PopulationInitializer::end() {
    return m_allele_freq.end();
}

std::vector< double >::const_iterator PopulationInitializer::end() const {
    return m_allele_freq.end();
}

PopulationInitializer::~PopulationInitializer() {}

}   // namespace genetics
}   // namespace clotho
