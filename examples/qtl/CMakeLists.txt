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

SET(LOCAL_LIBS)

LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_utility)
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_commandline)
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_fitness)
LIST(APPEND LOCAL_LIBS ${PROJECT_NAME}_popgrowth)
LIST(APPEND LOCAL_LIBS ${Boost_PROGRAM_OPTIONS_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_FILESYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_SYSTEM_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_CHRONO_LIBRARY})
LIST(APPEND LOCAL_LIBS ${COMMON_LIBRARIES} )

ADD_EXECUTABLE( qtl_estimator qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_64;USE_INLINE_DYNAMIC_AND_CLASSIFY" LINK_FLAGS "-Wl,--no-as-needed" )   # binary sequence combination determined by experimentation
TARGET_LINK_LIBRARIES( qtl_estimator ${LOCAL_LIBS} )

ADD_EXECUTABLE( qtl_estimator_vector qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator_vector PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_64;USE_VECTOR_SUBSET" LINK_FLAGS "-Wl,--no-as-needed" )
TARGET_LINK_LIBRARIES( qtl_estimator_vector ${LOCAL_LIBS})

#LIST(APPEND PROG_NAME "iterate")
#LIST(APPEND PROG_NAME "scan")
#LIST(APPEND PROG_NAME "switch")
#LIST(APPEND PROG_NAME "inline")
#LIST(APPEND PROG_NAME "inline_dynamic")
#LIST(APPEND PROG_NAME "vector")
#
#LIST(APPEND WALKER_OPT "USE_ITERATE_CLASSIFY")  # USE_ITERATE_CLASSIFY is never used with qtl compilation
#LIST(APPEND WALKER_OPT "USE_SCAN_AND_CLASSIFY")
#LIST(APPEND WALKER_OPT "USE_SCAN_AND_CLASSIFY_SWITCH")
#LIST(APPEND WALKER_OPT "USE_INLINE_AND_CLASSIFY")
#LIST(APPEND WALKER_OPT "USE_INLINE_DYNAMIC_AND_CLASSIFY")
#LIST(APPEND WALKER_OPT "USE_VECTOR_SUBSET")
#
#set(WORD_SIZE_OPT "USE_BLOCK_SIZE_")    # USE_BLOCK_SIZE_64 is never used in qtl compilation
#
#LIST(APPEND WORD_SIZE "32")
#LIST(APPEND WORD_SIZE "64")
#
#foreach( wsize ${WORD_SIZE} )
#    MESSAGE( STATUS "Got word size ${wsize}" )
#    SET( opt_idx 0)
#    foreach( pname  ${PROG_NAME} )
#        SET( tmp "qtl_estimator_${pname}_${wsize}")
#        LIST( GET WALKER_OPT ${opt_idx} tmp_walker)
#        SET( tmp_opts "${WORD_SIZE_OPT}${wsize};${tmp_walker}")
#        MESSAGE( STATUS "Generating Make for ${tmp} -> ${tmp_opts}")
#
#        ADD_EXECUTABLE( ${tmp} qtl_estimator.cpp )
#        SET_TARGET_PROPERTIES( ${tmp} PROPERTIES COMPILE_DEFINITIONS "${tmp_opts}")
#        TARGET_LINK_LIBRARIES( ${tmp} ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )
#        MATH(EXPR opt_idx "${opt_idx}+1")
#
#        LIST( APPEND TEST_TARGETS ${tmp} )
#    endforeach(pname)
#endforeach(wsize)
#
#ADD_CUSTOM_TARGET(test_targets DEPENDS ${TEST_TARGETS})

ADD_EXECUTABLE( qtl_estimator_mutrec qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator_mutrec PROPERTIES COMPILE_DEFINITIONS "USE_MUTATE_AND_RECOMBINE" LINK_FLAGS "-Wl,--no-as-needed")
TARGET_LINK_LIBRARIES( qtl_estimator_mutrec ${LOCAL_LIBS})

ADD_EXECUTABLE( qtl_estimator_assort qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator_assort PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_64;USE_VECTOR_SUBSET;USE_ASSORTATIVE_SELECTOR" LINK_FLAGS "-Wl,--no-as-needed")
TARGET_LINK_LIBRARIES( qtl_estimator_assort ${LOCAL_LIBS})
