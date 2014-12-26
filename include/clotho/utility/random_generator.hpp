#ifndef RANDOM_GENERATOR_HPP_
#define RANDOM_GENERATOR_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace utility {

template < class URNG, class Element >
class random_generator {
public:
    typedef random_generator< URNG, Element > self_type;
    typedef Element                           result_type;

    //random_generator( URNG & rng, boost::property_tree::ptree & config ) : m_rng( &rng ) {}

    result_type operator()() {
        return result_type( );
    }

protected:
    URNG * m_rng;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // RANDOM_GENERATOR_HPP_
