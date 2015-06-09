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
#ifndef DOUBLE_BUFFERED_MUTEX_STREAM_HPP_
#define DOUBLE_BUFFERED_MUTEX_STREAM_HPP_

#include "buffered_stream_def.hpp"
#include "clotho/utility/result_type_of.hpp"

#include <type_traits>

#include <pthread.h>
#include <iostream>

struct double_buffered_mutex {};

struct thread_data {
    pthread_mutex_t lock;
    pthread_cond_t  cond;

    bool            done;
    igenerator *    gen;
};

void * generator_worker_mutex( void * args ) {

    thread_data * data = (thread_data *)args;
    pthread_mutex_lock( &data->lock );

    while( true ) {
        pthread_cond_wait( &data->cond, &data->lock );

        if( data->done ) { break; }

        data->gen->generate();
    }

    pthread_mutex_unlock( &data->lock );

    pthread_exit(NULL);
}

template < class Generator >
class BufferedStream < Generator, double_buffered_mutex, typename std::enable_if< std::is_base_of< igenerator, Generator >::value >::type  > {
public:

    typedef typename clotho::utility::result_type_of< Generator >::type result_t;

    BufferedStream( Generator & gen, unsigned int size = 256) :
        m_buffer( NULL )
        , m_buffersize( size )
        , m_gen( &gen )
        , m_thread_err(false)
    {
        init();
    }

    result_t operator()() {
        if( m_cur == m_end ) {
            updateBuffer();
        }

        return *m_cur++;
    }

    bool threading_error() const {
        return m_thread_err;
    }
    
    virtual ~BufferedStream() {
        pthread_mutex_lock( &m_data.lock );

        m_data.done = true;

        pthread_cond_signal( &m_data.cond );
        pthread_mutex_unlock( &m_data.lock );

        pthread_join( m_worker, NULL );

        pthread_mutex_destroy( &m_data.lock );
        pthread_cond_destroy( &m_data.cond );
        pthread_attr_destroy( &m_attr );

        delete [] m_buffer;
    }
protected:
    void init() {
        // setup double buffer
        m_buffer = new result_t[ m_buffersize * 2 ];

        pthread_attr_init( &m_attr );
        pthread_mutex_init( &m_data.lock, NULL );
        pthread_cond_init( &m_data.cond, NULL );

        m_data.done = false;
        m_data.gen = (igenerator * )m_gen;

        // create worker thread
        int rc;
        rc = pthread_create( &m_worker, &m_attr, generator_worker_mutex, (void *) &m_data );
        m_thread_err = (rc != 0 );

        // fill initial buffer on main thread
        m_gen->start = m_buffer;
        m_gen->end = m_buffer + m_buffersize;

        m_gen->generate();

        // make it seem like we have reached the end
        // of the second buffer
        m_cur = m_buffer + 2 * m_buffersize;
        m_end = m_cur;
    }

    void updateBuffer( ) {

        pthread_mutex_lock( &m_data.lock );
        // fill
        if( m_end - m_buffer > m_buffersize ) {
            // switch to window 1
            m_cur = m_buffer;
            m_end = m_cur + m_buffersize;

            // [m_buffersize, 2 * m_buffersize)
            m_gen->start = m_end;
            m_gen->end = m_end + m_buffersize;
        } else {
            // continue to window 2
            m_end = m_cur + m_buffersize;

            // [ 0, m_buffersize)
            m_gen->start = m_buffer;
            m_gen->end = m_cur;
        }
        pthread_cond_signal( &m_data.cond );
        pthread_mutex_unlock( &m_data.lock );
    }

    result_t * m_buffer, *m_cur, *m_end;
    size_t m_buffersize;
    Generator * m_gen;
    bool            m_thread_err;
    thread_data     m_data;

    pthread_t       m_worker;
    pthread_attr_t  m_attr;
};

#endif // DOUBLE_BUFFERED_MUTEX_STREAM_HPP_
