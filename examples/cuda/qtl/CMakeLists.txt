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

LIST(APPEND LOCAL_LIBS ${CUDA_cuda_LIBRARY})
LIST(APPEND LOCAL_LIBS ${CUDA_curand_LIBRARY})
LIST(APPEND LOCAL_LIBS ${Boost_CHRONO_LIBRARY})
LIST(APPEND LOCAL_LIBS ${COMMON_LIBRARIES})

CUDA_ADD_EXECUTABLE(qtl_cuda qtl_cuda.cpp device_properties.h device_properties.cu  page_manager.cu adjacency_matrix.cu)
TARGET_LINK_LIBRARIES(qtl_cuda ${LOCAL_LIBS} )
