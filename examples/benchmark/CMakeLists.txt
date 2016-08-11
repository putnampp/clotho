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

#LIST( APPEND SRCS common_commandline.cpp )
#LIST( APPEND SRCS simulation_config_json.cpp )
#LIST( APPEND SRCS simulation_config.cpp )
#
#ADD_LIBRARY( ${PROJECT_NAME}_commandline SHARED ${SRCS} )

LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_utility)
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_commandline)
LIST(APPEND LOCAL_LIBS ${Boost_PROGRAM_OPTIONS_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_FILESYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_SYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_CHRONO_LIBRARY})

ADD_EXECUTABLE( benchmark_all_neutral benchmark.cpp )
SET_TARGET_PROPERTIES(benchmark_all_neutral PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES(benchmark_all_neutral ${LOCAL_LIBS})

ADD_EXECUTABLE( benchmark_selected benchmark.cpp )
SET_TARGET_PROPERTIES(benchmark_all_selected PROPERTIES COMPILE_DEFINITIONS "ALL_SELECTED_ALLELES" LINK_FLAGS "-Wl,--no-as-needed")
TARGET_LINK_LIBRARIES(benchmark_selected ${LOCAL_LIBS})
