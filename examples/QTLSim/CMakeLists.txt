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
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR}/qtl )

SET(LOCAL_LIBS)

LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_utility )
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_popgrowth )
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_config_manager )
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_fitness )

LIST(APPEND LOCAL_LIBS ${Boost_PROGRAM_OPTIONS_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_FILESYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_SYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_CHRONO_LIBRARY})

LIST(APPEND OPTIONS_SRC options/configuration_option.cpp )
LIST(APPEND OPTIONS_SRC options/log_prefix_option.cpp)
LIST(APPEND OPTIONS_SRC ../generation_parameter.cpp)
LIST(APPEND OPTIONS_SRC ../population_parameter.cpp)
LIST(APPEND OPTIONS_SRC ../logging_parameter.cpp)
LIST(APPEND OPTIONS_SRC ../qtl/qtl_logging_parameter.cpp)

LIST(APPEND OPTIONS_SRC simulation_log.cpp)
LIST(APPEND OPTIONS_SRC qtlsim_common.cpp)

IF( WITH_CUDA_EXAMPLES )

LIST(APPEND CUDA_LOCAL_LIBS ${LOCAL_LIBS} )
LIST(APPEND CUDA_LOCAL_LIBS ${CUDA_curand_LIBRARY})
LIST(APPEND CUDA_LOCAL_LIBS ${CUDA_cuda_LIBRARY})

LIST(APPEND CUDA_LOCAL_LIBS ${PROJECT_NAME}_cuda )
LIST(APPEND CUDA_LOCAL_LIBS ${PROJECT_NAME}_cuda_fitness )

SET(SIM_SRC main.cu)

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered ${SIM_SRC} ${OPTIONS_SRC} )
SET_TARGET_PROPERTIES( qtl_cuda_unordered PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered ${CUDA_LOCAL_LIBS})


CUDA_ADD_EXECUTABLE(qtl_cuda_unit_ordered ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DUSE_UNIT_ORDERING)
SET_TARGET_PROPERTIES( qtl_cuda_unit_ordered PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unit_ordered ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_verbose ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DDUMP_HEAP_VARIABLES)
SET_TARGET_PROPERTIES( qtl_cuda_unordered_verbose PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_verbose ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_double ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DUSE_DOUBLE_REAL)
SET_TARGET_PROPERTIES( qtl_cuda_unordered_double PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_double ${CUDA_LOCAL_LIBS})

ADD_CUSTOM_TARGET( qtl_cuda DEPENDS qtl_cuda_unordered qtl_cuda_unordered_verbose qtl_cuda_unordered_double )

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_offcross ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DUSE_OFFSPRING_CROSSOVER )
SET_TARGET_PROPERTIES( qtl_cuda_unordered_offcross PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_offcross ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_ham3 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DHAMMING_VERSION=3)
SET_TARGET_PROPERTIES( qtl_cuda_unordered_ham3 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_ham3 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_ham ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DHAMMING_VERSION)
SET_TARGET_PROPERTIES( qtl_cuda_unordered_ham PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_ham ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_ham4 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DHAMMING_VERSION=4)
SET_TARGET_PROPERTIES( qtl_cuda_unordered_ham4 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
#SET_TARGET_PROPERTIES( qtl_thrust PROPERTIES CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_ham4 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_xover2 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_VERSION=2 )
SET_TARGET_PROPERTIES( qtl_cuda_unordered_xover2 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_xover2 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_xover3 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_VERSION=3 )
SET_TARGET_PROPERTIES( qtl_cuda_unordered_xover3 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_xover3 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unordered_xover4 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_VERSION=4 )
SET_TARGET_PROPERTIES( qtl_cuda_unordered_xover4 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(qtl_cuda_unordered_xover4 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(qtl_cuda_unit_ordered_xover2 ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DUSE_UNIT_ORDERING;-DCROSSOVER_VERSION=2)
SET_TARGET_PROPERTIES( qtl_cuda_unit_ordered_xover2 PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(qtl_cuda_unit_ordered_xover2 ${CUDA_LOCAL_LIBS})

CUDA_ADD_EXECUTABLE(QTLSimCuda ${SIM_SRC} ${OPTIONS_SRC} OPTIONS -DCROSSOVER_VERSION=4)
SET_TARGET_PROPERTIES(QTLSimCuda PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(QTLSimCuda ${CUDA_LOCAL_LIBS})

ENDIF(WITH_CUDA_EXAMPLES)

LIST(APPEND LOCAL_LIBS ${COMMON_LIBRARIES})

SET(SIM_SRC main.cpp)

ADD_EXECUTABLE( QTLSim ${SIM_SRC} ${OPTIONS_SRC} )
SET_TARGET_PROPERTIES(QTLSim PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_64;USE_INLINE_DYNAMIC_AND_CLASSIFY;USE_FAST_GENERATOR" LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(QTLSim ${LOCAL_LIBS})