#ifndef IFITNESS_HPP_
#define IFITNESS_HPP_

#include <string>
#include <vector>
#include <ostream>

struct ifitness {
    typedef double result_type;

    /**
     * Reduce K-traits to a single 1D fitness value
     */
    virtual double operator()( const std::vector< double > & pheno ) = 0;
    virtual const std::string name() const = 0;

    virtual void log( std::ostream & out ) const = 0;

    virtual ~ifitness() {}
};

#endif  // IFITNESS_HPP_
