#   Copyright 2015 Patrick Putnam
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR}/cuda )

ADD_LIBRARY( cmdline_parser SHARED config_manager.cpp help_commandline_option.cpp)

SET(LOCAL_LIBS)

LIST(APPEND LOCAL_LIBS cmdline_parser )
LIST(APPEND LOCAL_LIBS ${CUDA_curand_LIBRARY})
LIST(APPEND LOCAL_LIBS ${CUDA_cuda_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_PROGRAM_OPTIONS_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_FILESYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_SYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_CHRONO_LIBRARY})
LIST(APPEND LOCAL_LIBS ${COMMON_LIBRARIES})

LIST(APPEND OPTIONS_SRC qtl_thrust_general_options.cpp)
LIST(APPEND OPTIONS_SRC random_number_options.cpp)
LIST(APPEND OPTIONS_SRC mutation_rate_option.cpp)
LIST(APPEND OPTIONS_SRC recombination_rate_option.cpp)
LIST(APPEND OPTIONS_SRC population_size_option.cpp)
LIST(APPEND OPTIONS_SRC generations_option.cpp)

LIST(APPEND SIM_SRC main.cpp)
LIST(APPEND SIM_SRC simulate_engine.cu )
LIST(APPEND SIM_SRC crossover_matrix.cu )
LIST(APPEND SIM_SRC population_recombiner.cu )
LIST(APPEND SIM_SRC scatter_mutations.cu )
LIST(APPEND SIM_SRC population_metadata.cu )

CUDA_ADD_EXECUTABLE(qtl_thrust ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_MATRIX=3)
SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_thrust ${LOCAL_LIBS} )

CUDA_ADD_EXECUTABLE(qtl_thrust_cm2 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_MATRIX=2)
SET_TARGET_PROPERTIES( qtl_thrust_cm2 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
TARGET_LINK_LIBRARIES(qtl_thrust_cm2 ${LOCAL_LIBS} )

# FEATURE TEST PROGRAM
CUDA_ADD_EXECUTABLE(kvsort_test kvsort_test.cu )
TARGET_LINK_LIBRARIES(kvsort_test ${CUDA_curand_LIBRARY} ${CUDA_cuda_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
