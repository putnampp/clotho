INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_SRC_DIR}/examples )

ADD_EXECUTABLE( qtl_estimator qtl_estimator.cpp )
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


#ADD_EXECUTABLE( qtl_estimator32 qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator32 PROPERTIES COMPILE_DEFINITIONS "USE_BLOCK_SIZE_32")
#TARGET_LINK_LIBRARIES( qtl_estimator32 ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )
#
#ADD_EXECUTABLE( qtl_estimator_bern qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_bern PROPERTIES COMPILE_DEFINITIONS "USE_BERNOULLI_RECOMB")
#TARGET_LINK_LIBRARIES( qtl_estimator_bern ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)
#
#ADD_EXECUTABLE( qtl_estimator_vector qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_vector PROPERTIES COMPILE_DEFINITIONS "USE_VECTOR_SUBSET")
#TARGET_LINK_LIBRARIES( qtl_estimator_vector ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)
#
#ADD_EXECUTABLE( qtl_estimator_scan_classify qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_scan_classify PROPERTIES COMPILE_DEFINITIONS "USE_SCAN_AND_CLASSIFY")
#TARGET_LINK_LIBRARIES( qtl_estimator_scan_classify ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)
#
#ADD_EXECUTABLE( qtl_estimator_scan_classify_switch qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_scan_classify_switch PROPERTIES COMPILE_DEFINITIONS "USE_SCAN_AND_CLASSIFY_SWITCH")
#TARGET_LINK_LIBRARIES( qtl_estimator_scan_classify_switch ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)
#
#ADD_EXECUTABLE( qtl_estimator_inline_classify qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_inline_classify PROPERTIES COMPILE_DEFINITIONS "USE_INLINE_AND_CLASSIFY")
#TARGET_LINK_LIBRARIES( qtl_estimator_inline_classify ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)
#
#ADD_EXECUTABLE( qtl_estimator_inline_dynamic_classify qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_inline_dynamic_classify PROPERTIES COMPILE_DEFINITIONS "USE_INLINE_DYNAMIC_AND_CLASSIFY")
#TARGET_LINK_LIBRARIES( qtl_estimator_inline_dynamic_classify ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric)

#
#ADD_EXECUTABLE( qtl_estimator_no_fit qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_no_fit PROPERTIES COMPILE_DEFINITIONS "CONSTFIT")
#TARGET_LINK_LIBRARIES( qtl_estimator_no_fit ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline )
