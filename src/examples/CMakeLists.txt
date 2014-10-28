INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_SRC_DIR}/examples )

LIST( APPEND SRCS common_commandline.cpp )
LIST( APPEND SRCS simulation_config_json.cpp )
LIST( APPEND SRCS simulation_config.cpp )

ADD_EXECUTABLE( benchmark_all_neutral benchmark.cpp ${SRCS} )
SET_TARGET_PROPERTIES(benchmark_all_neutral PROPERTIES COMPILE_DEFINITIONS "ALL_NEUTRAL_OPTIMIZATION")
TARGET_LINK_LIBRARIES(benchmark_all_neutral ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility)

ADD_EXECUTABLE( benchmark_selected benchmark.cpp ${SRCS} )
TARGET_LINK_LIBRARIES(benchmark_selected ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility)
