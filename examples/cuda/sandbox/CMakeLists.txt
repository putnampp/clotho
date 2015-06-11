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

ADD_EXECUTABLE( gen_test generate_test.cpp)

ADD_EXECUTABLE( buffer_test buffer_test.cpp)
TARGET_LINK_LIBRARIES( buffer_test ${COMMON_LIBRARIES} )
