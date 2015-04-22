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
