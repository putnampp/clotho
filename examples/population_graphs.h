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
#ifndef POPULATION_GRAPHS_H_
#define POPULATION_GRAPHS_H_

#include <vector>
#include <map>
#include <boost/property_tree/ptree.hpp>

#include "genetics/locus_bitset.h"

void logMutationFrequencies( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, std::map< size_t, size_t > & freq_dist, typename locus_bitset::alphabet_t::pointer alpha );

void logMutationBinFrequencies( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, unsigned int bin_count, typename locus_bitset::alphabet_t::pointer alpha );

void logMutationDistribution( boost::property_tree::ptree & p, std::map< size_t, size_t > & freq_dist );

void logMutationStats( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, typename locus_bitset::alphabet_t::pointer alpha );

#endif  // POPULATION_GRAPHS_H_
