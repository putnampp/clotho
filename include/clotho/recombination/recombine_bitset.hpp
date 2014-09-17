#ifndef RECOMBINE_BITSET_HPP_
#define RECOMBINE_BITSET_HPP_

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <limits>
#include <sstream>

#ifdef LOGGING
#include <boost/property_tree/ptree.hpp>
extern boost::property_tree::ptree global_log;
#endif

#include "utility/lowest_bit.h"

#include "accessor.hpp"

template < class Block, class Allocator, class AlleleSet >
class recombine_bitset;

namespace boost {

template < typename Block, typename Allocator, typename AlleleSet >
inline void to_block_range( const dynamic_bitset< Block, Allocator > & alt, recombine_bitset< Block, Allocator, AlleleSet > recomb ) {
    typedef typename std::vector< Block, Allocator >::const_iterator iterator;
    typedef Block block_type;

#ifdef LOGGING
    static unsigned int nCalls = 0;
    std::ostringstream oss;
    if( recomb.m_stats ) {
        oss << recomb.m_stats->gen << ".";
    }
    oss << "recombine_bitset." << nCalls++;
    std::string log_key = oss.str();

    global_log.put( log_key + ".alt.size", alt.size() );
    oss.str("");
    oss.clear();
    oss << alt;
    global_log.put( log_key + ".alt.sequence", oss.str());

    global_log.put( log_key + ".base.size", recomb.m_base->size());
    oss.str("");
    oss.clear();
    oss << *recomb.m_base;
    global_log.put( log_key + ".base.sequence", oss.str() );

    boost::property_tree::ptree pts;
    for( unsigned int i = 0; i < recomb.rec_points->size(); ++i ) {
        boost::property_tree::ptree t;
        t.put( "", recomb.rec_points->at(i));
        pts.push_back( std::make_pair( "", t ) );
    }

    global_log.add_child(log_key + ".points", pts );
#endif

    iterator first = recomb.m_base->m_bits.begin(), last = recomb.m_base->m_bits.end(),
             _first = alt.m_bits.begin(), _last = alt.m_bits.end();

    recomb( first, last, _first, _last );

    // truncate trailing bits
    recomb.m_result->resize( (alt.size() >= recomb.m_base->size()) ? alt.size() : recomb.m_base->size());

#ifdef LOGGING
    global_log.put( log_key + ".result.size", recomb.m_result->size());
    oss.str("");
    oss.clear();
    oss << *recomb.m_result;
    global_log.put( log_key + ".result.sequence", oss.str() );
    if( recomb.m_stats ) {
        global_log.put( log_key + ".match_base", recomb.m_stats->match_base );
        global_log.put( log_key + ".match_alt", recomb.m_stats->match_alt );
        global_log.put( log_key + ".is_empty", recomb.m_stats->is_empty );
    }

    if(recomb.m_base->count() ) {
        boost::property_tree::ptree a;
        for( typename dynamic_bitset< Block, Allocator >::size_type i = recomb.m_base->find_first(); i != dynamic_bitset< Block, Allocator >::npos; i = recomb.m_base->find_next(i)) {

            boost::property_tree::ptree t;
            t.put("", (*recomb.m_alpha)[i]->first);

            boost::property_tree::ptree r;
            typename recombine_bitset< Block, Allocator, AlleleSet >::recombination_iterator rit = std::upper_bound( recomb.rec_points->begin(), recomb.rec_points->end(), (*recomb.m_alpha)[i]->first);
            if( (rit - recomb.rec_points->begin()) % 2 ) {
                r.put("", "base");
            } else {
                r.put("", "alt");
            }

            boost::property_tree::ptree c;
            c.push_back( std::make_pair( "", t) );
            c.push_back( std::make_pair( "", r) );
            oss.str("");
            oss.clear();
            oss << i;
            a.add_child( oss.str(), c);
        }
        global_log.add_child( log_key + ".base.alpha", a);
    }

    if( alt.count() ) {
        boost::property_tree::ptree a;
        for( typename dynamic_bitset< Block, Allocator >::size_type i = alt.find_first(); i != dynamic_bitset< Block, Allocator >::npos; i = alt.find_next(i)) {

            boost::property_tree::ptree t;
            t.put("", (*recomb.m_alpha)[i]->first);

            boost::property_tree::ptree r;
            typename recombine_bitset< Block, Allocator, AlleleSet >::recombination_iterator rit = std::upper_bound( recomb.rec_points->begin(), recomb.rec_points->end(), (*recomb.m_alpha)[i]->first);
            if( (rit - recomb.rec_points->begin()) % 2 ) {
                r.put("", "base");
            } else {
                r.put("", "alt");
            }

            boost::property_tree::ptree c;
            c.push_back( std::make_pair( "", t) );
            c.push_back( std::make_pair( "", r) );
            oss.str("");
            oss.clear();
            oss << i;
            a.add_child( oss.str(), c);
        }
        global_log.add_child( log_key + ".alt.alpha", a);
    }
#endif

}

}   // namespace boost

template < class Block, class Allocator, class AlleleSet >
class recombine_bitset {
public:
    typedef recombine_bitset< Block, Allocator, AlleleSet > self_type;

    typedef boost::dynamic_bitset< Block, Allocator >   bitset_type;
    typedef std::vector< typename AlleleSet::locus_t >   recombination_points;
    typedef typename recombination_points::iterator     recombination_iterator;
    typedef typename AlleleSet::active_iterator          active_iterator;

    typedef lowest_bit_64K                              lowest_bit_map;

    struct result_stats {
        bool match_base, match_alt, is_empty;
#if LOGGING
        unsigned int gen;
        result_stats(unsigned int g = 0) : match_base(false), match_alt(false), is_empty(false), gen(g) {}
#else
        result_stats() : match_base(false), match_alt(false), is_empty(false) {}
#endif
    };

    typedef result_stats result_type;

    friend void boost::to_block_range< Block, Allocator, AlleleSet >( const boost::dynamic_bitset< Block, Allocator > &, recombine_bitset< Block, Allocator, AlleleSet>  );

    recombine_bitset( bitset_type * base, bitset_type * res,
                      typename AlleleSet::pointer alpha,
                      recombination_points * rp, result_type * stats = NULL ) :
        m_base( base )
        , m_result( res )
        , m_alpha( alpha )
        , rec_points( rp )
        , m_stats( stats ) {
        rec_points->insert(rec_points->begin(), std::numeric_limits< typename AlleleSet::locus_t >::min() );
        rec_points->push_back(std::numeric_limits< typename AlleleSet::locus_t >::max());
    }

    recombine_bitset( const self_type & other ) :
        m_base( other.m_base )
        , m_result( other.m_result )
        , m_alpha( other.m_alpha )
        , rec_points( other.rec_points )
        , m_stats( other.m_stats ) {
        rec_points->insert(rec_points->begin(), std::numeric_limits< typename AlleleSet::locus_t >::min() );
        rec_points->push_back(std::numeric_limits< typename AlleleSet::locus_t >::max());
    }

    template < class BlockIterator >
    void operator()( BlockIterator base_first, BlockIterator base_last, BlockIterator alt_first, BlockIterator alt_last ) {
        method1( base_first, base_last, alt_first, alt_last );
    }

    template < class BlockIterator >
    void method2( BlockIterator base_first, BlockIterator base_last, BlockIterator alt_first, BlockIterator alt_last ) {
        bool match_base = true, match_alt = true, is_empty = true;

        static bitset_type mask;
        mask.reset();
        mask.resize(m_alpha->size(), false);

        buildBaseMask( mask );

        BlockIterator mask_it = mask.m_bits.begin();

//        unsigned int nMaskBlocks = mask.num_blocks();

        while( true ) {
            if( base_first == base_last ) {
                while(alt_first != alt_last ) {
//                    assert(nMaskBlocks--);
                    Block a = (*alt_first++), m = (*mask_it++);

                    Block res = (a & ~m);

                    is_empty = ((is_empty) && (res == 0));
                    match_base = ((match_base) && (res == 0 ));
                    match_alt = ((match_alt) && (res == a));

                    m_result->append(res);
                }
                break;
            }

            if( alt_first == alt_last ) {
                while( base_first != base_last ) {
//                    assert(nMaskBlocks--);
                    Block b = (*base_first++), m = (*mask_it++);

                    Block res = (b & m);

                    is_empty = ((is_empty) && (res == 0));
                    match_base = ((match_base) && (res == b ));
                    match_alt = ((match_alt) && (res == 0 ));

                    m_result->append(res);
                }
                break;
            }

//            assert( nMaskBlocks-- );
            Block b = (*base_first++), a = (*alt_first++), m = (*mask_it++);

            Block res = ((b & m) | (a & ~m));

            is_empty = ((is_empty) && (res == 0));
            match_base = ((match_base) && (res == b ));
            match_alt = ((match_alt) && (res == a));

            m_result->append(res);
        }
    }

    template < class BlockIterator >
    void method1( BlockIterator base_first, BlockIterator base_last, BlockIterator alt_first, BlockIterator alt_last ) {
        bool match_base = true, match_alt = true, is_empty = true;

        active_iterator seq_pos = m_alpha->active_begin();
//        unsigned int seq_pos = 0;
//
//
        unsigned int nBase = (base_last - base_first), nAlt = (alt_last - alt_first);
        m_result->resize( (nBase > nAlt) ? (bitset_type::bits_per_block * nBase) : (bitset_type::bits_per_block * nAlt), false);

        typename std::vector< Block, Allocator>::iterator res_it = m_result->m_bits.begin();

        while( true ) {
            if( base_first == base_last ) {
                while( alt_first != alt_last ) {
                    Block alt = (*alt_first++);
                    Block res = alt;
                    bit_walker(res, alt, seq_pos, &self_type::unset_if_base);

                    is_empty = ((is_empty) && (res == 0));
                    match_base = ((match_base) && (res == 0 ));
                    match_alt = ((match_alt) && (res == alt));

//                    m_result->append(res);
                    (*res_it) = res;
                    ++res_it;
                    seq_pos += bitset_type::bits_per_block;
                }
                break;
            }

            if( alt_first == alt_last ) {
                while( base_first != base_last ) {
                    Block base = (*base_first++);
                    Block res = base;
                    bit_walker( res, base, seq_pos, &self_type::unset_if_alt );

                    is_empty = ((is_empty) && (res == 0));
                    match_base = ((match_base) && ( res == base ));
                    match_alt = ((match_alt) && (res == 0));

//                    m_result->append(res);
                    (*res_it) = res;
                    ++res_it;
                    seq_pos += bitset_type::bits_per_block;
                }
                break;
            }

            Block base = (*base_first++), alt = (*alt_first++);

            Block res = (base | alt);     // homozygous bits

            bit_walker( res, (base & ~alt), seq_pos, &self_type::unset_if_alt );
            bit_walker( res, (alt & ~base), seq_pos, &self_type::unset_if_base );

            is_empty = ((is_empty) && (res == 0));
            match_base = ((match_base) && (res == base ));
            match_alt = ((match_alt) && (res == alt ));

//            m_result->append( res );
            (*res_it) = res;
            ++res_it;
            seq_pos += bitset_type::bits_per_block;
        }

        if( m_stats ) {
            m_stats->match_base = match_base;
            m_stats->match_alt = match_alt;
            m_stats->is_empty = is_empty;
        }
    }

    virtual ~recombine_bitset() {}
protected:

    void buildBaseMask( bitset_type & mask ) {
        for( unsigned int i = 1; i < rec_points->size(); i += 2 ) {
            typename AlleleSet::locus_t min = rec_points->at(i - 1);
            typename AlleleSet::locus_t max = rec_points->at(i);

            typename AlleleSet::symbol_range res = m_alpha->findAllSymbols( min, max );

            while( res.first != res.second ) {
                if( res.first->second.second != bitset_type::npos )
                    mask[ res.first->second.second ] = true;
                ++res.first;
            }
        }
    }

    inline Block walk_block_bits( Block base, Block alt, unsigned int pos_offset ) {
        Block bits = (base ^ alt);      // symmetric difference
        Block res = (base & ~bits);     // homozygous bits

        Block bit_mask = 1;
        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            // NOTE: It may be more efficient to subdivide a block into byte regions
            // to allow for skipping of long regions of 0. Basically, recursively
            // walk the set bits by regions
            if( bits & 1 ) {    // if the locus is heterozygous
                // determine from which recombination region this locus is inherited from (base or alternate)
                recombination_iterator rit = std::lower_bound( rec_points->begin(), rec_points->end(), (*m_alpha)[pos_offset]->first);
                if( (rit - rec_points->begin()) % 2 ) {
                    res = ((res & ~bit_mask) | (alt & bit_mask));
                } else {
                    res = (( res & ~bit_mask) | (base & bit_mask));
                }
            }
            bits >>= 1;
            bit_mask <<= 1;
            ++pos_offset;
        }

        return res;
    }

    inline Block walk_block_bits_sparse( Block base, Block alt, unsigned int pos_offset ) {
        // performance should degrade as span of 0's is less than or equal to 8
        Block bits = (base ^ alt);      // symmetric difference
        Block res = (base & ~bits);     // homozygous bits

        Block bit_mask = 1;
        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            if( bits & 255 ) {
                unsigned int n = 8;
                while( n-- ) {
                    if( bits & 1 ) {    // if the locus is heterozygous
                        // determine from which recombination region this locus is inherited from (base or alternate)
                        recombination_iterator rit = std::lower_bound( rec_points->begin(), rec_points->end(), (*m_alpha)[pos_offset]->first);
                        if( (rit - rec_points->begin()) % 2 ) {
                            res = ((res & ~bit_mask) | (alt & bit_mask));
                        } else {
                            res = (( res & ~bit_mask) | (base & bit_mask));
                        }
                        n = 0;  // break out of loop to recheck state
                    }
                    bits >>= 1;
                    bit_mask <<= 1;
                    ++pos_offset;
                }
            } else {
                bits >>= 8;
                bit_mask <<= 8;
                pos_offset += 8;
            }
        }

        return res;
    }

    inline unsigned long walk_block_bits_sparse2( unsigned long base, unsigned long alt, unsigned int pos_offset ) {
        // performance should degrade as span of 0's is less than or equal to 8
        Block bits = (base ^ alt);      // symmetric difference
        Block res = (base & ~bits);     // homozygous bits

        Block bit_mask = 1;
        while( bits ) {     // foreach bit (locus)
            // this loop takes advantage of suffixes which are all 0
            // and truncates the computation after the last heterozygous position
            if( bits & 0x00000000000000FF ) {
                unsigned int n = 8;
                while( n-- ) {
                    if( bits & 1 ) {    // if the locus is heterozygous
                        // determine from which recombination region this locus is inherited from (base or alternate)
                        recombination_iterator rit = std::lower_bound( rec_points->begin(), rec_points->end(), (*m_alpha)[pos_offset]->first);
                        if( (rit - rec_points->begin()) % 2 ) {
                            res = ((res & ~bit_mask) | (alt & bit_mask));
                        } else {
                            res = (( res & ~bit_mask) | (base & bit_mask));
                        }
                        n = 0;  // break out of loop to recheck state
                    }
                    bits >>= 1;
                    bit_mask <<= 1;
                    ++pos_offset;
                }
            } else if( bits & 0x000000000000FF00) {
                bits >>= 8;
                bit_mask <<= 8;
                pos_offset += 8;
            } else if( bits & 0x00000000FFFF0000) {
                bits >>= 16;
                bit_mask <<= 16;
                pos_offset += 16;
            } else {
                bits >>= 32;
                bit_mask <<=32;
                pos_offset += 32;
            }
        }

        return res;
    }

    void unset_if_alt( Block & res, double loc, unsigned int _offset ) {
        recombination_iterator rit = std::upper_bound( rec_points->begin(), rec_points->end(), loc);
        if( (rit - rec_points->begin()) % 2 == 0 ) {
            // bit is supposed to come from alt
            // hence clear the bit in the result

#ifdef LOGGING
            static unsigned int n = 0;
            Block before = res;
#endif
            res &= ~(((Block)1) << (_offset));

#ifdef LOGGING
            if( n++ % 10 == 0) {
                std::ostringstream oss;
                oss << "unset_if_alt." << (n / 10);
                std::string log_key = oss.str();

                global_log.add( log_key + ".location", loc );
                global_log.add( log_key + ".upper_bound_offset", (rit - rec_points->begin()));
                global_log.add( log_key + ".upper_bound", *rit);

                oss.str("");
                oss.clear();
                oss << std::hex << (((Block)1) << (_offset));
                global_log.add( log_key + ".mask", oss.str());

                oss.str("");
                oss.clear();
                oss << std::hex << before;

                global_log.add( log_key + ".before", oss.str() );

                oss.str("");
                oss.clear();
                oss << std::hex << res;

                global_log.add( log_key + ".after", oss.str() );
            }
#endif
//        } else {
//            res |= (((Block)1) << (_offset));
        }
    }

    void unset_if_base( Block & res, double loc, unsigned int _offset ) {
//        double loc = (*m_alpha)[ lookup_offset + _offset ]->first;
        recombination_iterator rit = std::upper_bound( rec_points->begin(), rec_points->end(), loc);
        if( (rit - rec_points->begin()) % 2 ) {
            // bit is supposed to come from base
            // hence clear the bit in the result
#ifdef LOGGING
            static unsigned int n = 0;
            Block before = res;
#endif
            res &= ~(((Block)1) << (_offset));
#ifdef LOGGING
            if( n++ % 10 == 0) {
                std::ostringstream oss;
                oss << "unset_if_base." << (n / 10);
                std::string log_key = oss.str();

                global_log.add( log_key + ".location", loc );
                global_log.add( log_key + ".upper_bound_offset", (rit - rec_points->begin()));

                global_log.add( log_key + ".upper_bound", *rit);

                oss.str("");
                oss.clear();
                oss << std::hex << (((Block)1) << (_offset));
                global_log.add( log_key + ".mask", oss.str());

                oss.str("");
                oss.clear();
                oss << std::hex << before;

                global_log.add( log_key + ".before", oss.str() );

                oss.str("");
                oss.clear();
                oss << std::hex << res;

                global_log.add( log_key + ".after", oss.str() );
            }
#endif
//        } else {
//            res |= (((Block)1) << (_offset));
        }
    }

    typedef void (self_type::*unset_op_ptr)( Block & res, double, unsigned int boffset);

    /*    inline void bit_walker( Block & res, Block bits, unsigned int pos_offset,  unset_op_ptr op ) {
            unsigned int res_offset = 0;
            while( bits ) {
                unsigned char low_byte = (unsigned char)(bits & 0x00000000000000FF);
                unsigned int _offset = res_offset;
                while( low_byte ) {
                    const lowest_bit_map::value_type & v = low_bit_map[low_byte];
                    unsigned int shift = _offset + v.bit_index;

                    (this->*op)( res, (*m_alpha)[ pos_offset + shift]->first, shift );

                    low_byte = v.next;
                    _offset += v.bit_shift_next;
                }

                bits >>= 8;
                res_offset += 8;
            }
        }*/

    inline void rec_bit_walker( Block & res, unsigned short bits, active_iterator pos_offset, unset_op_ptr op, unsigned int res_offset ) {
        if( bits == 0 ) return;

        const lowest_bit_map::value_type * v = low_bit_map.begin() + bits;
        do {
            unsigned int shift = res_offset + v->bit_index;

//            (this->*op)( res, (*(pos_offset + shift))->first, shift);
            (this->*op)( res, accessor::get< typename AlleleSet::active_iterator, typename AlleleSet::locus_t >(pos_offset + shift), shift);

            res_offset += v->bit_shift_next;
            v = v->next_ptr;
        } while( v != NULL );
    }

    inline void rec_bit_walker( Block & res, unsigned int bits, active_iterator pos_offset, unset_op_ptr op, unsigned int res_offset ) {
        if( bits == 0 ) return;

        static const unsigned char WIDTH = 16;

        unsigned short _bits = (unsigned short) bits;
        rec_bit_walker( res, _bits, pos_offset, op, res_offset );

        _bits = (unsigned short)(bits >> WIDTH);

        rec_bit_walker( res, _bits, pos_offset, op, (res_offset + WIDTH) );
    }

    inline void rec_bit_walker( Block & res, unsigned long bits, active_iterator pos_offset, unset_op_ptr op ) {
        if( bits == 0 ) return;

        static const unsigned char WIDTH = 32;

        unsigned int _bits = (unsigned int) bits;
        rec_bit_walker( res, _bits, pos_offset, op, 0 );

        _bits = (unsigned int)(bits >> WIDTH);

        rec_bit_walker( res, _bits, pos_offset, op, WIDTH );
    }

    inline void bit_walker( Block & res, Block bits, active_iterator pos_offset,  unset_op_ptr op ) {
        unsigned int res_offset = 0;
        while( bits ) {
            typename lowest_bit_map::block_type low_byte = (typename lowest_bit_map::block_type)(bits);
            if( low_byte ) {
                unsigned int _offset = res_offset;
                const lowest_bit_map::value_type * v = low_bit_map.begin() + low_byte;
                do {
                    unsigned int shift = _offset + v->bit_index;

                    //(this->*op)( res, (*(pos_offset + shift))->first, shift);
                    (this->*op)( res, accessor::get< typename AlleleSet::active_iterator, typename AlleleSet::locus_t >( pos_offset + shift ), shift );

                    _offset += v->bit_shift_next;
                    v = v->next_ptr;
                } while( v != NULL );
            }

            bits >>= lowest_bit_map::block_width;
            res_offset += lowest_bit_map::block_width;
        }
    }

    inline Block walk_block_bits_sparse3( Block base, Block alt, unsigned int pos_offset ) {
        Block res = (base & ~( base ^ alt));     // homozygous bits

        bit_walker( res, (base & ~alt), pos_offset, &self_type::unset_if_alt );
        bit_walker( res, (alt & ~base), pos_offset, &self_type::unset_if_base );

        return res;
    }

    bitset_type * m_base, * m_result;
    typename AlleleSet::pointer  m_alpha;
    recombination_points * rec_points;
    result_type        * m_stats;

    static const lowest_bit_map low_bit_map;
};

template < class Block, class Allocator, class AlleleSet >
const typename recombine_bitset< Block, Allocator, AlleleSet >::lowest_bit_map recombine_bitset< Block, Allocator, AlleleSet >::low_bit_map;

#undef BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#endif  // RECOMBINE_BITSET_HPP_
