//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef FITNESS_EVALUATORS_HPP_
#define FITNESS_EVALUATORS_HPP_

/**
 *
 * Assume 1 compute block per individual
 *
 * Evaluates fitness from only first trait
 *
 */
template < class RealType >
__global__ void evaluate_quadratic_fitness( RealType * phenos, RealType * fitness, unsigned int trait_count, RealType stddev ) {
    assert( blockDimx.x * blockDim.y == 1 );

    // gridDim.x == individual count
    // 2 * gridDim.x == sequence count
    unsigned int p_idx = 2 * blockIdx.x;
    
    RealType fit = phenos[ p_idx ];
    fit += phenos[ p_idx + 1 ];

    fit /= stddev;
    fit *= fit;

    if( fit > 1.0) {
        fit = 0.0;
    } else {
        fit = 1.0 - fit;
    }

    fitness[ blockIdx.x ] = fit;
}


#endif  // FITNESS_EVALUATORS_HPP_
