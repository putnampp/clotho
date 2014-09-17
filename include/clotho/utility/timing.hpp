#ifndef TIMING_HPP_
#define TIMING_HPP_

#include <boost/chrono/chrono.hpp>
#include <boost/chrono/chrono_io.hpp>
#include <boost/chrono/system_clocks.hpp>

typedef boost::chrono::high_resolution_clock    clock_type;
typedef clock_type::time_point                  time_point_type;
typedef clock_type::duration                    duration_type;

/*!  Basic timer class */
class timer {
public:

    static constexpr double hertz = (double)clock_type::period::num/clock_type::period::den;
    timer() : m_start( clock_type::now() ), m_end( time_point_type::max() ) {}

    inline void start() {
        m_end = time_point_type::max();
        m_start = clock_type::now();
    }

    inline void stop() {
        m_end = clock_type::now();
    }

    duration_type elapsed() const {
        if( m_end == time_point_type::max() ) {
            time_point_type t = clock_type::now();
            return (t - m_start);
        }
        return (m_end - m_start);
    }

    virtual ~timer() {}
protected:
    time_point_type m_start, m_end;
};

#endif  // TIMING_HPP_
