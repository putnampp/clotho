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

struct genotype_distribution {
    unsigned int AA, AB, BA, BB;

    genotype_distribution( unsigned int aa, unsigned int ab, unsigned int ba, unsigned int bb ):
        AA(aa), AB(ab), BA(ba), BB(bb)
    {}

    genotype_distribution( const genotype_distribution & other ) :
        AA( other.AA ), AB( other.AB ), BA(other.BA), BB( other.BB )
    {}

    virtual ~genotype_distribution() {}
};

class PopulationInitializer {
public:

    typedef std::vector< genotype_distribution > genotype_dist_type;

    PopulationInitializer( boost::property_tree::ptree & config );

    bool hasDistribution() const;
    bool hasAlleleDistribution() const;
    bool hasGenotypeDistribution() const;

    size_t getAlleleCount() const;

    std::vector< double >::iterator allele_frequency_begin();
    std::vector< double >::const_iterator allele_frequency_begin() const;

    std::vector< double >::iterator allele_frequency_end();
    std::vector< double >::const_iterator allele_frequency_end() const;

    std::vector< genotype_distribution >::iterator genotype_begin();
    std::vector< genotype_distribution >::const_iterator genotype_begin() const;

    std::vector< genotype_distribution >::iterator genotype_end();
    std::vector< genotype_distribution >::const_iterator genotype_end() const;

    virtual ~PopulationInitializer();

private:
    std::vector< double > m_allele_freq;
    std::vector< genotype_distribution > m_geno_dist;
    unsigned int m_allele_count;
};

} // namespace genetics
} // namespace clotho

#endif  // CLOTHO_POPULATION_INITIALIZER_HPP_
