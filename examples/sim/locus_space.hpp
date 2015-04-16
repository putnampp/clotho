#ifndef LOCUS_SPACE_HPP_
#define LOCUS_SPACE_HPP_

#include "space.hpp"
#include "locus.hpp"

class locus_space : public space {
public:

    locus_space();

    virtual size_t size();
    virtual void prune();

    virtual ~locus_space();

protected:
    std::vector< std::shared_ptr< locus > > m_loci;
};

#endif  // LOCUS_SPACE_HPP_
