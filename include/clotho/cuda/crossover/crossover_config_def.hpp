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
#ifndef CROSSOVER_CONFIG_DEF_HPP_
#define CROSSOVER_CONFIG_DEF_HPP_

#include "clotho/cuda/crossover/xover_config.hpp"
#include "clotho/cuda/execution_configuration.hpp"
#include <typeinfo>

template < class OrderTag, unsigned char V >
struct kernel_exec< xover_config< OrderTag, V > > {
    typedef xover_config< OrderTag, V > xover_type;
    typedef kernel_exec< xover_type >   self_type;

    kernel_exec( boost::property_tree::ptree & config ) :
        _blocks( xover_type::MAX_BLOCKS )
        , _warps( xover_type::MAX_WARPS )
    {
        parse_configuration( config );
    }

    kernel_exec( const self_type & rhs ) :
        _blocks( rhs._blocks )
        , _warps( rhs._warps )
    {}

    void operator()( dim3 & bcount, dim3 & tcount ) {
        bcount.x = _blocks; // 1D grid
        bcount.y = 1;
        bcount.z = 1;

        tcount.x = 32;      // assume threads per warp
        tcount.y = _warps;
        tcount.z = 1;
    }

    void parse_configuration( boost::property_tree::ptree & config ) {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( "crossover.kernel", lconfig );

        _blocks = lconfig.get< unsigned int >( "blocks", xover_type::MAX_BLOCKS );
        _warps = lconfig.get< unsigned int >( "warps", xover_type::MAX_WARPS );

        _blocks = ((_blocks < xover_type::MAX_BLOCKS ) ? _blocks : xover_type::MAX_BLOCKS);
        _warps = ((_warps < xover_type::MAX_WARPS ) ? _warps : xover_type::MAX_WARPS );

        lconfig.put( "blocks", _blocks );
        lconfig.put( "warps", _warps );

        lconfig.put( "description.ordering", typeid( typename xover_type::order_tag_type ).name() );
        lconfig.put( "description.comment", "full warps only");
        lconfig.put( "description.algorithm_version", xover_type::VERSION );
        lconfig.put( "description.max_blocks", xover_type::MAX_BLOCKS );
        lconfig.put( "description.max_warps", xover_type::MAX_WARPS );

        config.put_child( "crossover.kernel", lconfig );
    }

    unsigned int _blocks, _warps;
};

#endif  // CROSSOVER_CONFIG_DEF_HPP_
