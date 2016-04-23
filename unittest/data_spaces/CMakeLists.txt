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

LIST(APPEND DATA_SPACE_SRCS allele_test.cpp )
LIST(APPEND DATA_SPACE_SRCS association_matrix_test.cpp )
LIST(APPEND DATA_SPACE_SRCS genetic_space_test.cpp )
LIST(APPEND DATA_SPACE_SRCS trait_accumulator_test.cpp )
LIST(APPEND DATA_SPACE_SRCS free_space_test.cpp )
LIST(APPEND DATA_SPACE_SRCS allele_generator_test.cpp )
LIST(APPEND DATA_SPACE_SRCS crossover_test.cpp )
LIST(APPEND DATA_SPACE_SRCS selection_test.cpp )
LIST(APPEND DATA_SPACE_SRCS mutation_test.cpp )


ADD_LIBRARY(${PROJECT_NAME}_data_spaces_test SHARED ${DATA_SPACE_SRCS} )
