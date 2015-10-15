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
#ifndef LOGGING_PARAMETER_HPP_
#define LOGGING_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>

extern const std::string LOG_BLOCK_K;
extern const std::string PERIOD_K;

struct logging_parameter {
    
    static const unsigned int DEFAULT_PERIOD = (unsigned int) -1;
    unsigned int m_period;

    logging_parameter( unsigned int p = DEFAULT_PERIOD ) :
        m_period( p ) 
    {}

    logging_parameter( boost::property_tree::ptree & config ) :
        m_period( DEFAULT_PERIOD )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( LOG_BLOCK_K, lconfig );

        m_period = lconfig.get< unsigned int >( PERIOD_K, m_period );

        lconfig.put( PERIOD_K, m_period);
        config.put_child( LOG_BLOCK_K, lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( PERIOD_K, m_period );
        l.put_child( LOG_BLOCK_K, c );
    }

    virtual ~logging_parameter() {}
};

#endif  // LOGGING_PARAMETER_HPP_
