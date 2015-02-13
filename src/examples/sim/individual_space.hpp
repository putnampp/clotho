#ifndef INDIVIDUAL_SPACE_HPP_
#define INDIVIDUAL_SPACE_HPP_

#include "individual.hpp"

#include <boost/property_tree/ptree.hpp>

class IndividualSpace {
public:

    IndividualSpace( boost::property_tree::ptree & config );

    std::shared_ptr< individual > create();

    size_t size();
    void prune();

    virtual ~IndividualSpace();
protected:
    virtual void init( boost::property_tree::ptree & config );

    std::vector< std::shared_ptr< individual > > m_individuals;
    individual_param    m_ind_param;
};

#endif  // INDIVIDUAL_SPACE_HPP_
