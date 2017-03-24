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
#ifndef QTLSIM_ENGINES_HPP_
#define QTLSIM_ENGINES_HPP_

#include "engine_def.hpp"

#ifdef USE_MT
#ifdef USE_PARALLEL_PIPELINE
#include "engine_parallel_pipeline.hpp"
#define ENGINE_MODE parallel_pipeline
#elif USE_PARALLEL_PIPELINE_BATCH
#include "engine_parallel_pipeline_batched.hpp"
#define ENGINE_MODE parallel_pipeline_batched
#else   // USE_PARALLEL_PIPELINE
#include "engine_batched_task.hpp"
#define ENGINE_MODE batched_task
#endif  // USE_PARALLEL_PIPELINE
#else   // USE_MT
#include "qtlsim_engine.hpp"
#define ENGINE_MODE base_engine
#endif  // USE_MT
#include "engine_batched_task.hpp"

#endif  // QTLSIM_ENGINES_HPP_
