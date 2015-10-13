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
#ifndef XOVER_CONFIG_HPP_
#define XOVER_CONFIG_HPP_


template < class OrderTag, unsigned char V >
struct xover_config {
    typedef OrderTag    order_tag_type;
    static const unsigned char VERSION;
    static const unsigned int MAX_BLOCKS;
    static const unsigned int MAX_WARPS;
};

template < class O, unsigned char V >
const unsigned char xover_config< O, V >::VERSION = V;

template < class O, unsigned char V >
const unsigned int xover_config< O, V >::MAX_BLOCKS = 10;

template < class O, unsigned char V >
const unsigned int xover_config< O, V >::MAX_WARPS = 4;

#include "clotho/cuda/data_spaces/tags/unordered_tag.hpp"

template <>
const unsigned int xover_config< unordered_tag, 3 >::MAX_WARPS = 16;


#include "clotho/cuda/data_spaces/tags/unit_ordered_tag.hpp"

template < class IntType >
struct xover_config< unit_ordered_tag< IntType >, 2 > {
    typedef unit_ordered_tag< IntType > order_tag_type;
    static const unsigned char          VERSION = 2;
    static const unsigned int           MAX_BLOCKS = 10;
    static const unsigned int           MAX_WARPS = 8;
};

#endif  // XOVER_CONFIG_HPP_
