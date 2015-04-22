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
#ifndef SIMULATION_CONFIG_JSON_H_
#define SIMULATION_CONFIG_JSON_H_

#include <string>
using std::string;

extern const string   CONFIG_BLOCK_K;

extern const string   RNG_BLOCK_K;
extern const string   GEN_BLOCK_K;
extern const string   POP_BLOCK_K;
extern const string   MUT_BLOCK_K;
extern const string   REC_BLOCK_K;
extern const string   REGION_BLOCK_K;
extern const string   LOG_BLOCK_K;

extern const string   MODEL_K;
extern const string   SIZE_K;
extern const string   SEED_K;
extern const string   GENERATOR_K;
extern const string   REPETITION_K;
extern const string   RATE_PER_REGION_K;
extern const string   RATE_PER_BASE_K;
extern const string   BASE_PER_REGION_K;
extern const string   PERIOD_K;

// optimizations
extern const string   OPT_BLOCK_K;
extern const string   CHECK_SELECTED_K;

// other keys
extern const string   MAX_K;
extern const string   MEAN_K;
extern const string   SIGMA_K;

#endif  // SIMULATION_CONFIG_JSON_H_
