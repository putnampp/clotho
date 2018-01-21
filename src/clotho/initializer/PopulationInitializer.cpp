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

PopulationInitializer::PopulationInitializer( boost::property_tree::ptree & config ) : m_allele_count(0) {
    boost::property_tree::ptree  init_tree, all_dist, geno_dist;
    init_tree = config.get_child( "initialize", init_tree);

    all_dist = init_tree.get_child( "allele_distribution", all_dist );
    for (auto& item : all_dist ) {
        m_allele_freq.push_back( item.second.get_value<double>() );
    }

    geno_dist = init_tree.get_child( "genotype_distribution", geno_dist );
    for (auto& item : geno_dist ) {
        genotype_distribution tmp( 0, 0, 0, 0);
        tmp.AA = item.second.get<unsigned int >( "AA", tmp.AA );
        tmp.AB = item.second.get<unsigned int >( "AB", tmp.AB );
        tmp.BA = item.second.get<unsigned int >( "BA", tmp.BA );
        tmp.BB = item.second.get<unsigned int >( "BB", tmp.BB );
        m_geno_dist.push_back( tmp );
    }

    if( m_allele_freq.size() > m_geno_dist.size() ) {
        m_allele_count = m_allele_freq.size();
    } else {
        m_allele_count = m_geno_dist.size();
    }
}

bool PopulationInitializer::hasDistribution() const {
    return !m_allele_freq.empty() || !m_geno_dist.empty();
}

bool PopulationInitializer::hasAlleleDistribution() const {
    return !m_allele_freq.empty();
}

bool PopulationInitializer::hasGenotypeDistribution() const {
    return !m_geno_dist.empty();
}

size_t PopulationInitializer::getAlleleCount() const {
    return m_allele_count;
}

std::vector< double >::iterator PopulationInitializer::allele_frequency_begin() {
    return m_allele_freq.begin();
}

std::vector< double >::const_iterator PopulationInitializer::allele_frequency_begin() const {
    return m_allele_freq.begin();
}

std::vector< double >::iterator PopulationInitializer::allele_frequency_end() {
    return m_allele_freq.end();
}

std::vector< double >::const_iterator PopulationInitializer::allele_frequency_end() const {
    return m_allele_freq.end();
}

std::vector< genotype_distribution >::iterator PopulationInitializer::genotype_begin() {
    return m_geno_dist.begin();
}

std::vector< genotype_distribution >::const_iterator PopulationInitializer::genotype_begin() const {
    return m_geno_dist.begin();
}

std::vector< genotype_distribution >::iterator PopulationInitializer::genotype_end() {
    return m_geno_dist.end();
}

std::vector< genotype_distribution >::const_iterator PopulationInitializer::genotype_end() const {
    return m_geno_dist.end();
}

PopulationInitializer::~PopulationInitializer() {}

}   // namespace genetics
}   // namespace clotho
