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

LIST( APPEND METRIC_SRCS normal_fitness_metric.cpp )
LIST( APPEND METRIC_SRCS normal_fitness_generator.cpp )
LIST( APPEND METRIC_SRCS quadratic_fitness_metric.cpp )
LIST( APPEND METRIC_SRCS quadratic_fitness_generator.cpp )
LIST( APPEND METRIC_SRCS fitness_toolkit.cpp )

ADD_LIBRARY( ${PROJECT_NAME}_fitmetric SHARED ${METRIC_SRCS} )

LIST( APPEND GROWTH_SRCS linear_population_growth.cpp )
LIST( APPEND GROWTH_SRCS linear_population_growth_generator.cpp )
LIST( APPEND GROWTH_SRCS predefined_population_growth.cpp )
LIST( APPEND GROWTH_SRCS predefined_population_growth_generator.cpp )
LIST( APPEND GROWTH_SRCS population_growth_toolkit.cpp )

ADD_LIBRARY( ${PROJECT_NAME}_popgrowth SHARED ${GROWTH_SRCS} )
