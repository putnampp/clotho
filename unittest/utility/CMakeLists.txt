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

LIST(APPEND UTILITY_SRCS bit_walker_test.cpp )
LIST(APPEND UTILITY_SRCS debruijn_bit_walker_test.cpp )
LIST(APPEND UTILITY_SRCS popcount_test.cpp )

ADD_LIBRARY(${PROJECT_NAME}_utility_test SHARED ${UTILITY_SRCS} )
