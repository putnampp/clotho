#ifndef PARAMETER_SPACE_HPP_
#define PARAMETER_SPACE_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

template < class Configable >
struct parameter_space {
    static void build_parameters( boost::property_tree::ptree & params ) {}
};

}   // namespace utility
}   // namespace clotho

#endif  // PARAMETER_SPACE_HPP_
