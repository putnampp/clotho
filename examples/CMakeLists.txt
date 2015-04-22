INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR} )

LIST( APPEND SRCS common_commandline.cpp )
LIST( APPEND SRCS simulation_config_json.cpp )
LIST( APPEND SRCS simulation_config.cpp )

ADD_LIBRARY( ${PROJECT_NAME}_commandline SHARED ${SRCS} )

#ADD_EXECUTABLE( benchmark_all_neutral benchmark.cpp )
#SET_TARGET_PROPERTIES(benchmark_all_neutral PROPERTIES COMPILE_DEFINITIONS "ALL_NEUTRAL_OPTIMIZATION")
#TARGET_LINK_LIBRARIES(benchmark_all_neutral ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline)
#
#ADD_EXECUTABLE( benchmark_selected benchmark.cpp )
#TARGET_LINK_LIBRARIES(benchmark_selected ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility ${PROJECT_NAME}_commandline)

ADD_SUBDIRECTORY(benchmark)
ADD_SUBDIRECTORY(qtl)
ADD_SUBDIRECTORY(cuda)
