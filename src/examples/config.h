#ifndef EXAMPLES_CONFIG_H_
#define EXAMPLES_CONFIG_H_

#ifdef ALL_NEUTRAL_OPTIMIZATION 
#define CHECK_SELECTED
#endif  // ALL_NEUTRAL_OPTIMIZATION

#include <string>
using std::string;

static const unsigned int DEFAULT_GENERATIONS = 10000;
static const unsigned int DEFAULT_POPULATION_SIZE = 10000;

static const double DEFAULT_MUTATION_RATE = 0.0001;
static const double DEFAULT_RECOMB_RATE = 0.0001;

static const unsigned int DEFAULT_SEED = 0;

#endif  // EXAMPLES_CONFIG_H_
