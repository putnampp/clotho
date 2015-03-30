#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <boost/chrono/chrono.hpp>
#include <boost/chrono/chrono_io.hpp>
#include <boost/chrono/system_clocks.hpp>

namespace clotho {
namespace utility {

typedef boost::chrono::high_resolution_clock    clock_type;
typedef clock_type::time_point                  time_point_type;
typedef clock_type::duration                    duration_type;

/*!  Basic timer class */
class timer {
public:

    static constexpr double hertz = (double)clock_type::period::num/clock_type::period::den;
    timer() : m_start( clock_type::now() ), m_end( time_point_type::max() ) {}

    /**
     * (re)start the timer
     */
    inline void start() {
        m_end = time_point_type::max();
        m_start = clock_type::now();
    }

    /**
     *  stop the timer
     */
    inline void stop() {
        m_end = clock_type::now();
    }

    /**
     * Calculates the elapsed time
     *
     * If the timer has been stopped, then the elasped time is stop - start.
     * Otherwise, the elapsed time is the time between "now" and the last start.
     *
     * \return elasped time
     * \note This does not stop a running timer
     */
    duration_type elapsed() const {
        time_point_type t = ((m_end == time_point_type::max()) ? clock_type::now() : m_end);
        return (t - m_start);
    }

    unsigned long elapsed_long() const {
        time_point_type t = ((m_end == time_point_type::max()) ? clock_type::now() : m_end);
        return (t - m_start).count();
    }

    virtual ~timer() {}
protected:
    time_point_type m_start, m_end;
};

}   // namespace utility
}   // namespace clotho
#endif  // TIMER_HPP_
