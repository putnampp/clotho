#ifndef INDIVIDUAL_FITNESS_HPP_
#define INDIVIDUAL_FITNESS_HPP_

template < class Fitness >
class individual_fitness {
public:
    typedef Fitness                             fitness_type;
    typedef typename fitness_type::set_type     set_type;
    typedef typename fitness_type::result_type  result_type;

    set_type *  m_set;

    individual_fitness( set_type & s ) : m_set(&s) {}

    template < class Individual >
    result_type operator()( Individual & ind ) {
        m_N += 1.;
        return (result_type) 1;
    }

    result_type total_fitness() const {
        return (result_type) m_N;
    }
    result_type expected_fitness() const {
        return (result_type) 1;
    }

protected:
    double m_N;
};



#include "clotho/powerset/variable_subset_fitness.hpp"
#include "neutral_method.hpp"

template < class E, class B, class BM, class EK, class Result, class HetMethod, class AltHomMethod, class RefHomMethod >
class individual_fitness<
    clotho::fitness::fitness< clotho::powersets::variable_subset< E, B, BM, EK >,
    HetMethod, AltHomMethod, RefHomMethod, Result > > {
public:
    typedef individual_fitness<
    clotho::fitness::fitness< clotho::powersets::variable_subset< E, B, BM, EK >,
           HetMethod, AltHomMethod, RefHomMethod, Result > > self_type;

    typedef clotho::fitness::fitness< clotho::powersets::variable_subset< E, B, BM, EK >,
            HetMethod, AltHomMethod, RefHomMethod, Result > fitness_type;
    typedef Result                                      result_type;
    typedef clotho::powersets::variable_subset< E, B, BM, EK >  sequence_type;
    typedef typename  sequence_type::powerset_type              set_type;
    typedef typename set_type::bitset_type                      bitset_type;

    //typedef neutral_method< E >                                 neutral_type;

    individual_fitness( set_type & s, const fitness_type & fit ) :
        m_set(&s)
        , m_fit( fit )
        , m_selected_mask()
        , m_default_value(fit.default_value())
        , m_tot_fit((result_type)0)
        , m_N((result_type)0)
        , m_nSelected(0) {
        initialize();
    }

    individual_fitness( const self_type & other ) :
        m_set( other.m_set )
        , m_fit( other.m_fit )
        , m_selected_mask( other.m_selected_mask )
        , m_default_value(other.m_default_value )
        , m_tot_fit(other.m_tot_fit )
        , m_N(other.m_N )
        , m_nSelected( other.m_nSelected ) {
    }

    template < class Individual >
    result_type operator()( Individual & ind );

    result_type operator()( std::pair< std::shared_ptr< sequence_type >, std::shared_ptr< sequence_type > > & ind ) {
        m_N += (result_type) 1;
        if(! m_nSelected ) {
            m_tot_fit += m_default_value;
            return m_default_value;
        }
        result_type res = m_fit( ind.first, ind.second, m_selected_mask.m_bits.begin(), m_selected_mask.m_bits.end() );
        m_tot_fit += res;

        return res;
    }

    result_type total_fitness() const {
        return m_tot_fit;
    }
    result_type expected_fitness() const {
        return m_tot_fit / ((m_N != 0) ? m_N : (result_type)1);
    }
    typename bitset_type::size_type selected_count() const {
        return m_nSelected;
    }

protected:
    void initialize() {
        m_selected_mask.clear();
        m_selected_mask.append( m_set->free_begin(), m_set->free_end() );
        m_selected_mask.flip();

        typename bitset_type::size_type idx = m_selected_mask.find_first();
        typename bitset_type::size_type offset = 0;
        unsigned int nSelected = 0;
        typename set_type::cvariable_iterator vit = m_set->variable_begin();
        while( idx != bitset_type::npos ) {
            offset = idx - offset;
            vit += offset;
            if( neutral_method::test( *vit ) ) {
                m_selected_mask.reset(idx);
            } else {
                ++nSelected;
            }
            offset = idx;
            idx = m_selected_mask.find_next(idx);
        }

        m_nSelected = nSelected;
    }

    set_type    * m_set;
    fitness_type    m_fit;
    bitset_type     m_selected_mask;
//    neutral_type    m_neutral;
    result_type     m_default_value, m_tot_fit, m_N;
    typename bitset_type::size_type m_nSelected;
};

#endif  // INDIVIDUAL_FITNESS_HPP_
