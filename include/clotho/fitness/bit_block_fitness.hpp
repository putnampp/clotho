#ifndef CLOTHO_BIT_BLOCK_FITNESS_HPP_
#define CLOTHO_BIT_BLOCK_FITNESS_HPP_

#include <iostream>
#include "clotho/fitness/bit_block_genotyper.hpp"
#include "clotho/utility/bit_block_iterator.hpp"

#include "clotho/fitness/no_fit.hpp"

namespace clotho {
namespace fitness {

template < class HetFit, class AltHomFit, class RefHomFit, class Result = double >
class bit_block_fitness {
public:
    typedef Result result_type;

    bit_block_fitness() {}

    bit_block_fitness( const HetFit & hets, const AltHomFit & ahoms, const RefHomFit & rhoms ) :
        m_hets( hets )
        , m_ahoms(ahoms)
        , m_rhoms(rhoms) {
    }

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        result_type res = f;

        typedef bit_block_genotyper< Block >                    genotyper;

        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_all_heterozygous, &m_hets);
        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_alt_homozygous, &m_ahoms);
        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_ref_homozygous, &m_rhoms);

        return res;
    }

    virtual ~bit_block_fitness() {}
protected:

    template < class Block, class ElementIterator, class FitOp >
    void computeFitness(result_type & res, Block b, ElementIterator first, FitOp * op ) {
        typedef clotho::utility::bit_block_iterator< Block, clotho::utility::walk_iterator_tag >    iterator;

        iterator it( b );
        while( !it.done()  ) {
            (*op)( res, *(first + *it));
            ++it;
        }
    }

    template < class Block, class ElementIterator, class SetOp, class FitOp >
    inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, FitOp * op ) {
        Block bits = ( (*sop)(b0, b1) & keep );
        computeFitness( res, bits, first, op );
    }

    template < class Block, class ElementIterator, class SetOp >
    inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, no_fit * op ) { }

    HetFit      m_hets;
    AltHomFit   m_ahoms;
    RefHomFit   m_rhoms;
};

template < class HetFit, class HomFit, class Result >
class bit_block_fitness< HetFit, HomFit, HomFit, Result > {
public:
    typedef Result result_type;

    bit_block_fitness( ) {}

    bit_block_fitness( const HetFit & hets, const HomFit & homs ) :
        m_hets( hets )
        , m_homs(homs) {
    }

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        result_type res = f;

        typedef bit_block_genotyper< Block >                    genotyper;

        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_all_heterozygous, &m_hets);
        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_alt_homozygous, &m_homs);

        return res;
    }

    virtual ~bit_block_fitness() {}
protected:
    template < class Block, class ElementIterator, class FitOp >
    inline void computeFitness(result_type & res, Block b, ElementIterator first, FitOp * op ) {
        typedef clotho::utility::bit_block_iterator< Block >    iterator;

        iterator it( b );
        while( !it.done() ) {
            (*op)( res, *(first + *it));
            ++it;
        }
    }

    template < class Block, class ElementIterator, class SetOp, class FitOp >
    inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, FitOp * op ) {
        Block bits = ( (*sop)(b0, b1) & keep );
        computeFitness( res, bits, first, op );
    }

    template < class Block, class ElementIterator, class SetOp >
    inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, no_fit * op ) { }

    HetFit   m_hets;
    HomFit   m_homs;
};

template < class Result >
class bit_block_fitness< no_fit, no_fit, no_fit, Result > {
public:
    typedef Result result_type;

    bit_block_fitness( ) {}

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        return f;
    }

    virtual ~bit_block_fitness() {}
};

}   // namespace fitness {
}   // namespace clotho {

#endif  // CLOTHO_BIT_BLOCK_FITNESS_HPP_
