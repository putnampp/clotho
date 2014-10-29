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

extern const string   MODEL_K;
extern const string   SIZE_K;
extern const string   SEED_K;
extern const string   GENERATOR_K;
extern const string   REPETITION_K;
extern const string   RATE_PER_REGION_K;
extern const string   RATE_PER_BASE_K;
extern const string   BASE_PER_REGION_K;
extern const string   LOG_FREQUENCY_K;

// optimizations
extern const string   OPT_BLOCK_K;
extern const string   CHECK_SELECTED_K;

#endif  // SIMULATION_CONFIG_JSON_H_
