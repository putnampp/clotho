#ifndef ALLELE_HPP_
#define ALLELE_HPP_

#include <memory>

#include "locus.hpp"
#include 

class Allele {
public:
    typedef unsigned int    allele_state_type;

    Allele();

    std::shared_ptr< locus >        getLocus();
    std::shared_ptr< effect_size >  getEffectSize();

    virtual ~Allele();

protected:
    std::shared_ptr< locus >        m_locus;
    std::shared_ptr< effect_size >  m_effect;
    allele_state_type               m_state;
};

#endif  // ALLELE_HPP_
