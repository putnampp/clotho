#ifndef ALLELE_SPACE_HPP_
#define ALLELE_SPACE_HPP_

#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "allele.hpp"

/**
 * Allele Space maps implicitly to Effect Space.  That is, the allele
 * at index i of AlleleSpace maps to the effect_size at index i in
 * EffectSpace.
 *
 * Alleles maps explicitly to loci in LocusSpace.
 *
 */
class AlleleSpace {
public:
    AlleleSpace( boost::property_tree::ptree & config );

    std::shared_ptr< Allele > create();

    size_t size();
    void prune();

    virtual ~AlleleSpace();
protected:
    EffectSpace                                 m_effects;
    LocusSpace                                  m_loci;
    std::vector< std::shared_ptr< Allele > >    m_alleles;

};

#endif  // ALLELE_SPACE_HPP_
