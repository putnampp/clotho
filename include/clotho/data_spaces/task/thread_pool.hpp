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
template < class RNG >
class thread_pool {
public:

    typedef RNG     random_engine_type;

    thread_pool( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_available( 0 )
        , m_pool_size( 0 )
        , m_running( true )
        , m_thread_engines( NULL )
    {
        thread_count_parameter tc_param( config );
        m_pool_size = tc_param.m_tc;

        
        m_thread_engines = new random_engine_type[ m_pool_size + 1];

        for( size_t i = 0; i < m_pool_size; ++i ) {
            m_threads.create_thread( boost::bind( &thread_pool::thread_loop, this ) );

            // seed each engine with a random number from the parent random number generator
            unsigned long long seed = (*rng)();

            BOOST_LOG_TRIVIAL(info) << "Thread " << i << " seed value: " << seed;
            m_thread_engines[ i ] = random_engine_type( seed );
        }
    }

    random_engine_type * getRNG( unsigned int index ) {
#ifdef DEBUGGING
        assert( index < m_pool_size );
#endif  // DEBUGGING
        return &m_thread_engines[ index ];
    }

    size_t pool_size() const {
        return m_pool_size;
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

        delete [] m_thread_engines;
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

    // one random number generator per thread
    random_engine_type                      * m_thread_engines;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_THREAD_POOL_HPP_

