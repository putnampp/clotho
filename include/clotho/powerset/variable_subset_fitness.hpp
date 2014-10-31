#ifndef CLOTHO_VARIABLE_SUBSET_FITNESS_HPP_
#define CLOTHO_VARIABLE_SUBSET_FITNESS_HPP_

#include "clotho/fitness/fitness_def.hpp"
#include "clotho/fitness/bit_block_fitness.hpp"

namespace clotho {
namespace fitness {

template < class Element, class Block, class BlockMap, class ElementKeyer, class HetFit, class AltHomFit, class RefHomFit, class Result >
class fitness< clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer >, HetFit, AltHomFit, RefHomFit, Result > {
public:
    typedef clotho::powersets::variable_subset< Element, Block, BlockMap, ElementKeyer > subset_type;
    typedef typename subset_type::pointer   sequence_type;
    typedef Result  result_type;
    typedef Block   block_type;

    typedef clotho::fitness::bit_block_fitness< HetFit, AltHomFit, RefHomFit, Result >  fitness_type;

    static const unsigned int bits_per_block = sizeof(Block) * 8;
    
    result_type operator()( result_type f, sequence_type base, sequence_type alt ) {
        if( base == alt ) {
            // pointers match
            return f;   // null sequences
        }

        typename subset_type::cblock_iterator base_it, base_end;
        typename subset_type::cblock_iterator alt_it, alt_end;

        typename subset_type::powerset_type::cvariable_iterator elem_it;

        if( !base ) {
            // base sequence is empty
            // NOTE: since base != alt, alt must be defined
            alt_it = alt->begin();
            alt_end = alt->end();

            // slight "trick" to cause following loop to iterate only
            // over the alt sequence
            base_it = alt->end();
            base_end = alt->end();

            elem_it = alt->getParent()->variable_begin();
        } else if( !alt ) {
            // alt sequence is empty
            base_it = base->begin();
            base_end = base->end();

            // slight "trick" to cause following loop to iterate only
            // over the base sequence
            alt_it = base->end();
            alt_end = base->end();

            elem_it = base->getParent()->variable_begin();
        } else {
            assert( base->isSameFamily( alt ) );

            base_it = base->begin();
            base_end = base->end();

            alt_it = alt->begin();
            alt_end = alt->end();

            elem_it = base->getParent()->variable_begin();
        }

        fitness_type       bfitness;

        block_type mask = (block_type) -1;
        result_type res = f;
        while( true ) {
            if( alt_it == alt_end ) {
                while( base_it != base_end ) {
                    block_type _base = (*base_it++);
                    res = bfitness( res, _base, (block_type)0, mask, elem_it );
                    elem_it += bits_per_block;
                }
                break;
            }

            if( base_it == base_end ) {
                while( alt_it != alt_end ) {
                    block_type _alt = (*alt_it++);
                    res = bfitness(res, (block_type)0, _alt, mask, elem_it );
                    elem_it += bits_per_block;
                }
                break;
            }

            block_type _base = (*base_it++), _alt = (*alt_it++);
            res = bfitness(res, _base, _alt, mask, elem_it );
            elem_it += bits_per_block;
        }

        return res;
    }

};

}   // namespace fitness
}   // namespace clotho

#endif  // CLOTHO_VARIABLE_SUBSET_FITNESS_HPP_
