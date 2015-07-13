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
#ifndef CLOTHO_TIMER_HPP_
#define CLOTHO_TIMER_HPP_

#include <boost/chrono/chrono.hpp>
#include <boost/chrono/chrono_io.hpp>
#include <boost/chrono/system_clocks.hpp>
#include <ostream>

namespace clotho {
namespace utility {

typedef boost::chrono::high_resolution_clock    clock_type;
typedef clock_type::time_point                  time_point_type;
typedef clock_type::duration                    duration_type;

/*!  Basic timer class */
class timer {
public:

    //static constexpr double hertz = (double)clock_type::period::num/(double)clock_type::period::den;
    static double hertz() {
        static double h = (double) clock_type::period::num / (double) clock_type::period::den;
        return h;
    }

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

    friend std::ostream & operator<<( std::ostream & out, const timer & t );

    virtual ~timer() {}
protected:
    time_point_type m_start, m_end;
};

inline std::ostream & operator<<( std::ostream & out, const timer & t ) {
    out << t.elapsed().count();
    return out;
}

}   // namespace utility
}   // namespace clotho
#endif  // CLOTHO_TIMER_HPP_
