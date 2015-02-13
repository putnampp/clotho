#ifndef SEQUENCE_SPACE_HPP_
#define SEQUENCE_SPACE_HPP_

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/dynamic_bitset.hpp>

#include "allele_space.hpp"

class SequenceSpace {
public:
    typedef boost::dynamic_bitset< unsigned long >  sequence_type; 
    typedef std::shared_ptr< sequence_type >        sequence_ptr;

    SequenceSpace( boost::property_tree::ptree & config );

    sequence_ptr create();

    sequence_ptr clone( sequence_ptr s );

    size_t size();
    void prune();

    virtual ~SequenceSpace();
protected:
    AlleleSpace m_alleles;

    std::vector< sequence_ptr > m_seq_map;
};

#endif  // SEQUENCE_SPACE_HPP_
