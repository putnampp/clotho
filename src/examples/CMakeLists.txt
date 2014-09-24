INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_SRC_DIR}/examples )


#LIST( APPEND SRCS common_strings.cpp )
LIST( APPEND SRCS common_commandline.cpp )
LIST( APPEND SRCS population_graphs.cpp )

ADD_EXECUTABLE(benchmark_all_neutral benchmark_all_neutral.cpp ${SRCS} )
TARGET_LINK_LIBRARIES( benchmark_all_neutral ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility)

