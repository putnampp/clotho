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
#ifndef QTL_LOGGING_PARAMETER_HPP_
#define QTL_LOGGING_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>
#include "../common_strings.h"
#include "../logging_parameter.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

extern const string SAMPLING_K;
extern const string PAIRWISE_K;

struct sample_log_params {
    static const unsigned int   DEFAULT_SAMPLE_SIZE = 100;
    static const bool           DEFAULT_PAIRWISE = false;
    unsigned int    sample_size;
    bool            pairwise;

    sample_log_params( unsigned int ss = DEFAULT_SAMPLE_SIZE, bool p = DEFAULT_PAIRWISE ) :
        sample_size(ss)
        , pairwise(p) 
    {}

    sample_log_params( boost::property_tree::ptree & params ) :
        sample_size(DEFAULT_SAMPLE_SIZE)
        , pairwise(DEFAULT_PAIRWISE)
    {
        if( params.empty() ) {
            // done for backwards compatiability
            // Check for data in the node
            std::ostringstream tmp;
            tmp << params.data();

            if(! tmp.str().empty() ) {
                // assume that any data was intended to be a sample size
                sample_size = boost::lexical_cast< unsigned int >( tmp.str() );
            }
            return;
        }

        sample_size = params.get< unsigned int >(SIZE_K, sample_size);

        pairwise = params.get< bool >( PAIRWISE_K, pairwise);
        params.put( SIZE_K, sample_size);
        params.put( PAIRWISE_K, pairwise );
    }

    sample_log_params( const sample_log_params & slp ) :
        sample_size( slp.sample_size )
        , pairwise( slp.pairwise ) {
    }

    virtual ~sample_log_params() {}
};

struct qtl_logging_parameter : public logging_parameter {

    std::vector< sample_log_params > m_sampling;

    qtl_logging_parameter( boost::property_tree::ptree & config ) : 
        logging_parameter( config )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( LOG_BLOCK_K, lconfig );

        boost::property_tree::ptree samps;
        samps = lconfig.get_child( SAMPLING_K, samps );

        BOOST_FOREACH( auto& v, samps ) {
            sample_log_params tmp( v.second );
            m_sampling.push_back( tmp );
        }

        lconfig.put_child( SAMPLING_K, samps );
        config.put_child( LOG_BLOCK_K, lconfig );
    }

    virtual ~qtl_logging_parameter() {}
};
#endif  // QTL_LOGGING_PARAMETER_HPP_
