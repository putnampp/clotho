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
#ifndef THREAD_COUNT_PARAMETER_HPP_
#define THREAD_COUNT_PARAMETER_HPP_

#include <boost/property_tree/ptree.hpp>

struct thread_count_parameter {
    int    m_tc;

    static const int DEFAULT_THREAD_COUNT = 0;

    thread_count_parameter( int m = DEFAULT_THREAD_COUNT ) :
        m_tc( m )
    { }

    thread_count_parameter( boost::property_tree::ptree & config ) :
        m_tc( DEFAULT_THREAD_COUNT )
    {
        boost::property_tree::ptree lconfig;
        lconfig = config.get_child( "multithreading", lconfig );

        m_tc = lconfig.get< int >( "T", m_tc );

        lconfig.put( "T", m_tc );
        config.put_child( "multithreading", lconfig );
    }

    void write_parameter( boost::property_tree::ptree & l ) {
        boost::property_tree::ptree c;
        c.put( "worker_threads", m_tc );
        l.put_child( "multithreading", c );
    }

    virtual ~thread_count_parameter() {}
};

#endif  // THREAD_COUNT_PARAMETER_HPP_
