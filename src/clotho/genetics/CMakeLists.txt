INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_INCLUDE_DIR} )

LIST( APPEND METRIC_SRCS normal_fitness_metric.cpp )
LIST( APPEND METRIC_SRCS normal_fitness_generator.cpp )
LIST( APPEND METRIC_SRCS quadratic_fitness_metric.cpp )
LIST( APPEND METRIC_SRCS quadratic_fitness_generator.cpp )
LIST( APPEND METRIC_SRCS fitness_toolkit.cpp )

ADD_LIBRARY( ${PROJECT_NAME}_fitmetric SHARED ${METRIC_SRCS} )
