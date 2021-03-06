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
#ifndef CLOTHO_THREAD_POOL2_HPP_
#define CLOTHO_THREAD_POOL2_HPP_

#include <vector>
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
class thread_pool2 {
public:

    typedef RNG     random_engine_type;

    thread_pool2( random_engine_type * rng, boost::property_tree::ptree & config ) :
        m_available( 0 )
        , m_pool_size( 0 )
        , m_running( true )
        , m_thread_engines( NULL )
    {
        thread_count_parameter tc_param( config );
        m_pool_size = tc_param.m_tc;

        
        m_thread_engines = new random_engine_type[ m_pool_size + 1];

        for( size_t i = 0; i < m_pool_size; ++i ) {

            // seed each engine with a random number from the parent random number generator
            unsigned long long seed = (*rng)();

            BOOST_LOG_TRIVIAL(info) << "Thread " << i << " seed value: " << seed;
            m_thread_engines[ i ] = random_engine_type( seed );
        }
    }

    thread_pool2( unsigned int T ) :
        m_threads(NULL)
        , m_available(0)
        , m_pool_size( T )
        , m_running( true )
        , m_thread_engines( NULL ) 
    {
        m_threads = new boost::thread* [ m_pool_size ];
        for( unsigned int i = 0; i < m_pool_size; ++i ) {
            m_threads[ i ] = NULL;
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

//    template < class Task >
//    void post( Task & t ) {
//        boost::thread *th = new boost::thread( t );
//
//        m_threads.add_thread( th );
//    }

    template < class Task >
    void post_list( std::vector< Task * > & task ) {
        assert( task.size() <= m_pool_size );

        unsigned int i = 0;
        for(typename std::vector< Task * >::iterator it = task.begin(); it != task.end() && i < m_pool_size; it++ ) {
//            boost::thread * th = new boost::thread( &Task::operator(), (*it) );
//            m_threads.add_thread( th );
            assert( m_threads[ i ] == NULL );
            m_threads[ i++ ] = new boost::thread( &Task::operator(), (*it));           
        }
    }

    template < class Task >
    void post_list( Task ** t, const unsigned int N ) {
        assert( N <= m_pool_size );
        for( unsigned int i = 0; i < N; ++i ) {
//            boost::thread * th = new boost::thread( &Task::operator(), t[i] );
//            m_threads.add_thread( th );
            assert( m_threads[ i ] == NULL );
            m_threads[ i ] = new boost::thread( &Task::operator(), t[ i ] );
        }
    }

    void sync() {
//        m_threads.join_all();
        for( unsigned int i = 0; i < m_pool_size; ++i ) {
            if( m_threads[ i ] != NULL ) {
                boost::thread * th = m_threads[ i ];
                m_threads[ i ] = NULL;
                th->join();
                delete th;
            }
        }
    }

    virtual ~thread_pool2() {
//        m_threads.join_all();
        sync();

        delete [] m_threads;

        if( m_thread_engines != NULL ) {
            delete [] m_thread_engines;
        }
    }

protected:

    boost::thread   ** m_threads;  
//    boost::thread_group                     m_threads;
    size_t                                  m_available, m_pool_size;
    bool                                    m_running;

    // one random number generator per thread
    random_engine_type                      * m_thread_engines;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_THREAD_POOL2_HPP_

