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
__global__ void evaluate_quadratic_fitness( RealType * phenos, RealType * fitness, RealType stddev ) {
    assert( blockDim.x * blockDim.y == 1 );

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

/**
 * 
 * Assume 1 thread per individual
 * Phenotype are aligned by trait
 */
template < class RealType >
__global__ void evaluate_quadratic_fitness( RealType * phenos, RealType * fitness, unsigned int individual_count, RealType stddev ) {

    unsigned int p_idx = threadIdx.y * blockDim.x + threadIdx.x;

    while( p_idx < individual_count ) {
        RealType fit = phenos[ 2 * p_idx ];
        fit += phenos[ 2 * p_idx + 1 ];

        fit /= stddev;
        fit *= fit;

        if( fit > 1.0) {
            fit = 0.0;
        } else {
            fit = 1.0 - fit;
        }

        fitness[ p_idx ] = fit;
        p_idx += blockDim.x * blockDim.y;
    }
}

#endif  // FITNESS_EVALUATORS_HPP_
