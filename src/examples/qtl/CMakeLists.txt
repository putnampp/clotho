INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_SRC_DIR}/examples )

ADD_EXECUTABLE( qtl_estimator qtl_estimator.cpp )
TARGET_LINK_LIBRARIES( qtl_estimator ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline ${PROJECT_NAME}_fitmetric )

#ADD_EXECUTABLE( qtl_estimator_quad qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_quad PROPERTIES COMPILE_DEFINITIONS "QFIT")
#TARGET_LINK_LIBRARIES( qtl_estimator_quad ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline )
#
#ADD_EXECUTABLE( qtl_estimator_no_fit qtl_estimator.cpp )
#SET_TARGET_PROPERTIES(qtl_estimator_no_fit PROPERTIES COMPILE_DEFINITIONS "CONSTFIT")
#TARGET_LINK_LIBRARIES( qtl_estimator_no_fit ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline )
