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
#include "simulation_config_json.h"

//const string   CONFIG_BLOCK_K = "configuration";

//const string   RNG_BLOCK_K = "random_number";
//const string   GEN_BLOCK_K = "generations";
//const string   POP_BLOCK_K = "population";
//const string   MUT_BLOCK_K = "mutation";
//const string   REC_BLOCK_K = "recombination";
const string   REGION_BLOCK_K = "region";
//const string   LOG_BLOCK_K = "log";

const string   MODEL_K = "model";
//const string   SIZE_K = "size";
//const string   SEED_K = "seed";
//const string   REPETITION_K = "repeat";
const string   GENERATOR_K = "generator";
const string   RATE_PER_REGION_K = "rate_per_region";
const string   RATE_PER_BASE_K = "rate_per_base";
const string   BASE_PER_REGION_K = "base_per_region";
//const string   PERIOD_K = "period";

const string   OPT_BLOCK_K = "compiled_optimizations";
const string   CHECK_SELECTED_K = "check_selected";

const string   MAX_K = "max";
const string   MEAN_K = "mean";
const string   SIGMA_K = "sigma";
