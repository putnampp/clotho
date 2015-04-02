#include <boost/python.hpp>

using namespace boost::python;

#include "clotho/utility/timer.hpp"

BOOST_PYTHON_MODULE( my_timer ) {
    class_< clotho::utility::timer >( "timer" )
        .def( "start", &clotho::utility::timer::start )
        .def( "stop", &clotho::utility::timer::stop )
        .def( "elapsed_long", &clotho::utility::timer::elapsed_long)
    ;
}
