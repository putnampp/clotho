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

LIST( APPEND SRCS common_commandline.cpp )
LIST( APPEND SRCS simulation_config_json.cpp )
LIST( APPEND SRCS simulation_config.cpp )

ADD_LIBRARY( ${PROJECT_NAME}_commandline SHARED ${SRCS} )

#ADD_EXECUTABLE( benchmark_all_neutral benchmark.cpp )
#SET_TARGET_PROPERTIES(benchmark_all_neutral PROPERTIES COMPILE_DEFINITIONS "ALL_NEUTRAL_OPTIMIZATION")
#TARGET_LINK_LIBRARIES(benchmark_all_neutral ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline)
#
#ADD_EXECUTABLE( benchmark_selected benchmark.cpp )
#TARGET_LINK_LIBRARIES(benchmark_selected ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline)

ADD_SUBDIRECTORY(benchmark)
ADD_SUBDIRECTORY(qtl)
ADD_SUBDIRECTORY(cuda)
