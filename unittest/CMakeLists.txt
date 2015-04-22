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
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIR})

LIST(APPEND DIRS powerset )
LIST(APPEND DIRS utility )
LIST(APPEND DIRS fitness )

FOREACH( D ${DIRS} )
    ADD_SUBDIRECTORY( ${D} )
ENDFOREACH( D )

ADD_EXECUTABLE(test_main test_main.cpp)
TARGET_LINK_LIBRARIES( test_main ${COMMON_LIBRARIES} ${PROJECT_NAME}_powerset_test ${PROJECT_NAME}_utility ${PROJECT_NAME}_utility_test ${PROJECT_NAME}_fitness_test)
