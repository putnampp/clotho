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
