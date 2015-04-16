INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR}/examples )

ADD_EXECUTABLE( qtl_estimator qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_64;USE_INLINE_DYNAMIC_AND_CLASSIFY" )   # binary sequence combination determined by experimentation
TARGET_LINK_LIBRARIES( qtl_estimator ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )

LIST(APPEND PROG_NAME "iterate")
LIST(APPEND PROG_NAME "scan")
LIST(APPEND PROG_NAME "switch")
LIST(APPEND PROG_NAME "inline")
LIST(APPEND PROG_NAME "inline_dynamic")
LIST(APPEND PROG_NAME "vector")

LIST(APPEND WALKER_OPT "USE_ITERATE_CLASSIFY")  # USE_ITERATE_CLASSIFY is never used with qtl compilation
LIST(APPEND WALKER_OPT "USE_SCAN_AND_CLASSIFY")
LIST(APPEND WALKER_OPT "USE_SCAN_AND_CLASSIFY_SWITCH")
LIST(APPEND WALKER_OPT "USE_INLINE_AND_CLASSIFY")
LIST(APPEND WALKER_OPT "USE_INLINE_DYNAMIC_AND_CLASSIFY")
LIST(APPEND WALKER_OPT "USE_VECTOR_SUBSET")

set(WORD_SIZE_OPT "USE_BLOCK_SIZE_")    # USE_BLOCK_SIZE_64 is never used in qtl compilation

LIST(APPEND WORD_SIZE "32")
LIST(APPEND WORD_SIZE "64")

foreach( wsize ${WORD_SIZE} )
    MESSAGE( STATUS "Got word size ${wsize}" )
    SET( opt_idx 0)
    foreach( pname  ${PROG_NAME} )
        SET( tmp "qtl_estimator_${pname}_${wsize}")
        LIST( GET WALKER_OPT ${opt_idx} tmp_walker)
        SET( tmp_opts "${WORD_SIZE_OPT}${wsize};${tmp_walker}")
        MESSAGE( STATUS "Generating Make for ${tmp} -> ${tmp_opts}")

        ADD_EXECUTABLE( ${tmp} qtl_estimator.cpp )
        SET_TARGET_PROPERTIES( ${tmp} PROPERTIES COMPILE_DEFINITIONS "${tmp_opts}")
        TARGET_LINK_LIBRARIES( ${tmp} ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )
        MATH(EXPR opt_idx "${opt_idx}+1")

        LIST( APPEND TEST_TARGETS ${tmp} )
    endforeach(pname)
endforeach(wsize)

ADD_CUSTOM_TARGET(test_targets DEPENDS ${TEST_TARGETS})

ADD_EXECUTABLE( qtl_estimator_mutrec qtl_estimator.cpp )
SET_TARGET_PROPERTIES( qtl_estimator_mutrec PROPERTIES COMPILE_DEFINITIONS "USE_MUTATE_AND_RECOMBINE")
TARGET_LINK_LIBRARIES( qtl_estimator_mutrec ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )
