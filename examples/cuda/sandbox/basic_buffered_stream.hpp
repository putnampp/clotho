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
#ifndef BASIC_BUFFERED_STREAM_HPP_
#define BASIC_BUFFERED_STREAM_HPP_

#include "buffered_stream_def.hpp"
#include <type_traits>

#include "clotho/utility/result_type_of.hpp"

struct basic_buffer {};

template < class Generator >
class BufferedStream< Generator, basic_buffer, typename std::enable_if< std::is_base_of< igenerator, Generator >::value >::type > {
public:
    typedef typename clotho::utility::result_type_of< Generator >::type result_t;

    BufferedStream( Generator & gen, unsigned int size = 256 ) :
        m_buffer(NULL)
        , m_buffersize( size )
        , m_gen( & gen )
    {
        init();
    }

    result_t operator()() {
        if( m_cur == m_end ) {
            updateBuffer();
        }
        return *m_cur++;
    }

    virtual ~BufferedStream() {
        delete [] m_buffer;
    }
protected:

    void init() {
        // setup double buffer
        m_buffer = new result_t[ m_buffersize ];

        m_cur = m_buffer;
        m_end = m_cur + m_buffersize;

        // fill buffer
        fillBuffer( m_cur, m_end );
    }

    void updateBuffer( ) {
        fillBuffer( m_buffer, m_end );
        m_cur = m_buffer;
    }

    void fillBuffer( result_t * s, result_t * e ) {
        m_gen->start =(void *) s;
        m_gen->end = (void *) e;
        m_gen->generate();
    }

    result_t * m_buffer, *m_cur, *m_end;
    size_t m_buffersize;
    Generator * m_gen;
};

#endif  // BASIC_BUFFERED_STREAM_HPP_
