#ifndef STATE_OF_HPP_
#define STATE_OF_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

template < class Object >
struct state_of {

    static void record( Object & obj, boost::property_tree::ptree & s ) {}
};

}   // namespace utility
}   // namespace clotho

#endif  // STATE_OF_HPP_
