#ifndef CLOTHO_BIT_BLOCK_FITNESS2_HPP_
#define CLOTHO_BIT_BLOCK_FITNESS2_HPP_

#include <iostream>
#include "clotho/fitness/bit_block_genotyper.hpp"
#include "clotho/utility/bit_block_iterator.hpp"

#include "clotho/fitness/no_fit.hpp"

namespace clotho {
namespace fitness {

template < class HetFit, class AltHomFit, class RefHomFit, class Result = double >
class bit_block_fitness2 {
public:
    typedef Result result_type;

    bit_block_fitness2() {}

    bit_block_fitness2( const HetFit & hets, const AltHomFit & ahoms, const RefHomFit & rhoms ) :
        m_hets( hets )
        , m_ahoms(ahoms)
        , m_rhoms(rhoms) {
    }

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        result_type res = f;

        typedef bit_block_genotyper< Block >                    genotyper;

//        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_all_heterozygous, &m_hets);
//        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_alt_homozygous, &m_ahoms);
//        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_ref_homozygous, &m_rhoms);
//

        Block het_bits = genotyper::get_all_heterozygous(b0, b1) & keep_mask;
        Block hom_bits = genotyper::get_alt_homozygous(b0, b1) & keep_mask;
        Block ref_bits = genotyper::get_ref_homozygous(b0, b1) & keep_mask;

//        Block _bits = (het_bits | (hom_bits | ref_bits) );
        Block _bits = 0;

        combineBits( _bits, het_bits, &m_hets );
        combineBits( _bits, hom_bits, &m_ahoms );
        combineBits( _bits, ref_bits, &m_rhoms );

        typedef clotho::utility::bit_block_iterator< Block, clotho::utility::walk_iterator_tag >    iterator;
        iterator it( _bits );

        while( !it.done() ) {
            unsigned int idx = (*it);
            ++it;
            Block _m = ((Block)1 << idx);

            if( het_bits & _m ) {
                //m_hets( res, *(first + idx) );
                computeFitness( res, *(first + idx), &m_hets );
            } else if( hom_bits & _m ) {
                //m_ahoms( res, *(first + idx) );
                computeFitness( res, *(first + idx), &m_ahoms );
            } else {
                computeFitness( res, *(first + idx), &m_rhoms );
            }
        }

        return res;
    }

    virtual ~bit_block_fitness2() {}
protected:

    template < class Element, class FitOp >
    inline void computeFitness( result_type & res, const Element & elem, FitOp * op ) {
        (*op)(res, elem);
    }

    template < class Element >
    inline void computeFitness( result_type & res, const Element & elem, no_fit * op ) { }

    template < class Block, class FitOp >
    inline void combineBits( Block & res, Block & set, FitOp * op ) {
        res |= set;
    }

    template < class Block >
    inline void combineBits( Block & res, Block & set, no_fit * op ) {}

    /*
        template < class Block, class ElementIterator, class FitOp >
        void computeFitness(result_type & res, Block b, ElementIterator first, FitOp * op ) {
            typedef clotho::utility::bit_block_iterator< Block, clotho::utility::walk_iterator_tag >    iterator;

            iterator it( b ), end;
            while( it != end ) {
                (*op)( res, *(first + *it++));
            }
        }

        template < class Block, class ElementIterator, class SetOp, class FitOp >
        inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, FitOp * op ) {
            Block bits = ( (*sop)(b0, b1) & keep );
            computeFitness( res, bits, first, op );
        }

        template < class Block, class ElementIterator, class SetOp >
        inline void computeFitness( result_type & res, Block b0, Block b1, Block keep, ElementIterator first, SetOp * sop, no_fit * op ) { }
    */
    HetFit      m_hets;
    AltHomFit   m_ahoms;
    RefHomFit   m_rhoms;
};

template < class HetFit, class HomFit, class Result >
class bit_block_fitness2< HetFit, HomFit, HomFit, Result > {
public:
    typedef Result result_type;

    bit_block_fitness2( ) {}

    bit_block_fitness2( const HetFit & hets, const HomFit & homs ) :
        m_hets( hets )
        , m_homs(homs) {
    }

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        result_type res = f;

        typedef bit_block_genotyper< Block >                    genotyper;

        std::cerr << "Bad call" << std::endl;

        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_all_heterozygous, &m_hets);
        computeFitness( res, b0, b1, keep_mask, first, &genotyper::get_alt_homozygous, &m_homs);

        return res;
    }

    virtual ~bit_block_fitness2() {}
protected:
    template < class Block, class ElementIterator, class FitOp >
    inline void computeFitness(result_type & res, Block b, ElementIterator first, FitOp * op ) {
        typedef clotho::utility::bit_block_iterator< Block, clotho::utility::walk_iterator_tag >    iterator;

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
class bit_block_fitness2< no_fit, no_fit, no_fit, Result > {
public:
    typedef Result result_type;

    bit_block_fitness2( ) {}

    template < class Block, class ElementIterator >
    result_type operator()( result_type f, Block b0, Block b1, Block keep_mask, ElementIterator first ) {
        return f;
    }

    virtual ~bit_block_fitness2() {}
};

}   // namespace fitness2 {
}   // namespace clotho {

#endif  // CLOTHO_BIT_BLOCK_FITNESS2_HPP_
