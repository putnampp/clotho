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
LIST(APPEND DIRS classifiers )

SET(LOCAL_LIBS)

LIST(APPEND LOCAL_LIBS ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

FOREACH( D ${DIRS} )
    ADD_SUBDIRECTORY( ${D} )
    LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_${D}_test)
ENDFOREACH( D )

LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_utility )
LIST(APPEND LOCAL_LIBS ${COMMON_LIBRARIES})

MESSAGE( STATUS "Unittest library list: " ${LOCAL_LIBS})

ADD_EXECUTABLE(test_main test_main.cpp)
SET_TARGET_PROPERTIES( test_main PROPERTIES LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES( test_main ${LOCAL_LIBS})
