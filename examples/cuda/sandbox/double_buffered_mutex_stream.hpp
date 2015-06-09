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
#ifndef DOUBLE_BUFFERED_STREAM_HPP_
#define DOUBLE_BUFFERED_STREAM_HPP_

#include "buffered_stream_def.hpp"
#include "clotho/utility/result_type_of.hpp"

#include <type_traits>

#include <pthread.h>

struct double_buffered {};

void * generator_worker( void * args ) {

    ((igenerator *)args)->generate();

    pthread_exit(NULL);
}

template < class Generator >
class BufferedStream < Generator, double_buffered, typename std::enable_if< std::is_base_of< igenerator, Generator >::value >::type  > {
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
        pthread_join( m_worker, NULL );

        delete [] m_buffer;
    }
protected:
    void init() {
        // setup double buffer
        m_buffer = new result_t[ m_buffersize * 2 ];

        pthread_attr_init( &m_attr );

        // fill first buffer
        fillBuffer( m_buffer, m_buffer + m_buffersize );

        // make it seem like we have reached the end
        // of the second buffer
        m_cur = m_buffer + 2 * m_buffersize;
        m_end = m_cur;
    }

    void updateBuffer( ) {
        int rc;
        void *status;

        if( !m_thread_err ) {
            rc = pthread_join( m_worker, &status );
            m_thread_err = (rc != 0);
        }
        // fill
        if( m_end - m_buffer > m_buffersize ) {
            // switch to window 1
            m_cur = m_buffer;
            m_end = m_cur + m_buffersize;
            fillBuffer( m_end, m_end + m_buffersize );   // [m_buffersize, 2 * m_buffersize)
        } else {
            // continue to window 2
            m_end = m_cur + m_buffersize;
            fillBuffer( m_buffer, m_cur );  // [ 0, m_buffersize)
        }
    }

    void fillBuffer( result_t * s, result_t * e) {
        m_gen->start = s;
        m_gen->end = e;

        int rc;
        rc = pthread_create( &m_worker, &m_attr, generator_worker, (void *) m_gen );
        m_thread_err = (rc != 0 );
    }

    result_t * m_buffer, *m_cur, *m_end;
    size_t m_buffersize;
    Generator * m_gen;
    bool            m_thread_err;

    pthread_t       m_worker;
    pthread_attr_t  m_attr;
};

#endif // DOUBLE_BUFFERED_STREAM_HPP_
