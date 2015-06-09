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
#include <iostream>

#include "clotho/utility/timer.hpp"

#include "buffered_stream.hpp"

typedef clotho::utility::timer timer_type;

class simple_generator : public igenerator {
public:
    typedef unsigned int result_t;

    simple_generator( unsigned int o = 0 ) : igenerator( NULL, NULL ), offset(o) {}

    void generate() {
        if( start == NULL || end == NULL ) return;

        result_t * t = static_cast< result_t *>( start );
        result_t * ee = static_cast< result_t *>( end );

        while( t < ee ) {
            *t++ = offset++;
        }
    }

    virtual ~simple_generator() {}

protected:
    unsigned int offset;
};

namespace clotho {
namespace utility {

template <>
struct result_type_of< simple_generator > {
    typedef unsigned int type;
};

}   // namespace utility
}   // namespace clotho

template < class Func >
void runTimedFunc( Func & f ) {
    timer_type t;
    f();
    t.stop();
    std::cout << "Lapse: " << t.elapsed().count() << std::endl;
}

template < class Stream, class Func >
void runTest( Stream & s, Func f, unsigned int L ) {
    std::cout << "Before: " << f << std::endl;

    timer_type t;
    while( L-- ) {
        f( s() );
    }
    t.stop();

    std::cout << "Result: " << f << std::endl;
    std::cout << "Lapse: " << t.elapsed().count() << std::endl;
}

template < class T >
struct summation {
    T val;

    summation( T v = 0) : val( v ) {}
    summation( const summation< T > & rhs ) : val(rhs.val) {}

    void operator()( T s ) { val += s; }
};

template < class T >
std::ostream & operator<<( std::ostream & out, const summation< T > & rhs ) {
    out << rhs.val;
    return out;
}

template < class Stream >
struct stream_summation {
    typedef typename clotho::utility::result_type_of< Stream >::type result_t;
    Stream * m_stream;
    unsigned int L;

    stream_summation( Stream & s, unsigned int l ) : m_stream( &s ), L(l) {}

    void operator()() {
        result_t res = 0;
        while( L-- ) {
            res += (*m_stream)();   
        }
        std::cout << "Result: " << res << std::endl;
    }
};

int main( int argc, char ** argv ) {

    unsigned int N = 1000000;

    typedef simple_generator                                            Generator;
    typedef BufferedStream< Generator, basic_buffer >                   SimpleBuffer;
    typedef BufferedStream< Generator, double_buffered >                DoubleBuffer;

    {
    Generator s(0);
    SimpleBuffer simple_bs(s, N/4);

    stream_summation< SimpleBuffer > test( simple_bs, N );
//    runTest( simple_bs, summation< result_type >(), N );
        runTimedFunc( test );
    }

    {
    Generator s2(0);
    DoubleBuffer double_bs(s2, N/5);

//    runTest( double_bs, summation< result_type >(), N );
        stream_summation< DoubleBuffer > test( double_bs, N );
        runTimedFunc(test);
    }

    return 0;
}
