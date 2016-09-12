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
#ifndef CLOTHO_THREAD_POOL_HPP_
#define CLOTHO_THREAD_POOL_HPP_

#include <queue>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include <boost/property_tree/ptree.hpp>
#include "clotho/data_spaces/task/thread_count_parameter.hpp"

namespace clotho {
namespace genetics {

/**
 *
 * Initial object construction from: http://stackoverflow.com/questions/12215395/thread-pool-using-boost-asio
 */
class thread_pool {
public:

    thread_pool( boost::property_tree::ptree & config ) :
        m_available( 0 )
        , m_pool_size( 0 )
        , m_running( true )
    {
        thread_count_parameter tc_param( config );
        m_pool_size = tc_param.m_tc;

        for( size_t i = 0; i < m_pool_size; ++i ) {
            m_threads.create_thread( boost::bind( &thread_pool::thread_loop, this ) );
        }
    }

    template < class Task >
    void post( Task t ) {
        boost::unique_lock< boost::mutex > lock( m_mutex );

        m_tasks.push( boost::function< void() >( t ) );

        if( 0 < m_available ) {
            m_condition.notify_one();
        }
    }

    void sync() {
        boost::unique_lock< boost::mutex > lock( m_mutex );
        if( !m_tasks.empty() || m_available < m_pool_size ) {
            m_sync_cond.wait( lock );
        }
    }

    virtual ~thread_pool() {
        {
            boost::unique_lock< boost::mutex > lock( m_mutex );
            m_running = false;
            m_condition.notify_all();
        }

        m_threads.join_all();
    }

protected:

    void thread_loop() {
        while( m_running ) {
            boost::unique_lock< boost::mutex > lock( m_mutex );
            while( m_tasks.empty() && m_running ) {
                if( ++m_available == m_pool_size ) {
                    m_sync_cond.notify_one();
                }
                m_condition.wait( lock );
                --m_available;
            }

            if( !m_running ) break;

            {
                boost::function< void() > f = m_tasks.front();
                m_tasks.pop();
                lock.unlock();

                f();
            }
        }
    }

    std::queue< boost::function< void() > > m_tasks;
    boost::thread_group                     m_threads;
    size_t                                  m_available, m_pool_size;
    boost::mutex                            m_mutex;
    boost::condition_variable             m_condition, m_sync_cond;
    bool                                    m_running;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_THREAD_POOL_HPP_

