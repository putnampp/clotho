#ifndef TOOLKIT_HPP_
#define TOOLKIT_HPP_

#include "clotho/utility/igenerator.hpp"

#include <boost/property_tree/ptree.hpp>

template < class T >
class toolkit {
public:

    std::shared_ptr< T > create_generator( boost::property_tree::ptree & config );
};

#endif  // TOOLKIT_HPP_
