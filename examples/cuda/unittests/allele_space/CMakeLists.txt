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
INCLUDE_DIRECTORIES( ${${PROJECT_NAME}_EXAMPLE_DIR}/cuda )

# FEATURE TEST PROGRAM
#CUDA_ADD_EXECUTABLE(allele_space_unordered_test allele_space_test.cu)
#TARGET_LINK_LIBRARIES(allele_space_unordered_test ${CUDA_curand_LIBRARY} ${CUDA_cuda_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_CHRONO_LIBRARY})

#CUDA_ADD_EXECUTABLE(allele_space_unit_ordered_test allele_space_test.cu OPTIONS -DUSE_UNIT_ORDERED)
#TARGET_LINK_LIBRARIES(allele_space_unit_ordered_test ${CUDA_curand_LIBRARY} ${CUDA_cuda_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_CHRONO_LIBRARY})

#ADD_CUSTOM_TARGET( allele_space_tests DEPENDS allele_space_unit_ordered_test allele_space_unordered_test )
