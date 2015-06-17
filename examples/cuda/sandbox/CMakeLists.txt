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
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_SRC_DIR}/examples )


CUDA_ADD_EXECUTABLE(square_list square_list.cpp square.cu square_host.cu cuda_mt19937.cu square_stream.cu)
TARGET_LINK_LIBRARIES(square_list ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY} ${COMMON_LIBRARIES} ${PROJECT_NAME}_utility )
LIST( APPEND SANDBOX square_list )

ADD_EXECUTABLE( gen_test generate_test.cpp)
LIST( APPEND SANDBOX gen_test )

ADD_EXECUTABLE( buffer_test_simple buffer_test.cpp)
TARGET_LINK_LIBRARIES( buffer_test_simple ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY} ${COMMON_LIBRARIES} )
LIST( APPEND SANDBOX buffer_test_simple )

CUDA_ADD_EXECUTABLE( buffer_test_cuda buffer_test.cpp)
SET_TARGET_PROPERTIES( buffer_test_cuda PROPERTIES COMPILE_DEFINITIONS "USE_CUDA_GENERATOR")
TARGET_LINK_LIBRARIES( buffer_test_cuda ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY} ${COMMON_LIBRARIES} )
LIST( APPEND SANDBOX buffer_test_cuda )

ADD_CUSTOM_TARGET( sandbox DEPENDS ${SANDBOX} )
