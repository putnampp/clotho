#ifndef PAIRWISE_STATISTIC_VECTOR_HPP_
#define PAIRWISE_STATISTIC_VECTOR_HPP_

#include "clotho/genetics/pairwise_statistic_def.hpp"

#include "clotho/powerset/vector_subset.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/weighted_sum.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

namespace accum=boost::accumulators;

template < class E, class B, class BM, class EK >
struct pairwise_statistic< clotho::powersets::vector_subset< E, B, BM, EK > > {

    typedef clotho::powersets::vector_subset< E, B, BM, EK > sequence_type;
    typedef accum::accumulator_set< double, accum::stats< accum::tag::weighted_mean, accum::tag::weighted_variance, accum::tag::weighted_sum, accum::tag::count >, double > accum_type;

    pairwise_statistic( sequence_type & b, double w = 1.0 ) :
        _base( &b )
        , base_weight(w)
        , m_global_diff( NULL )
        , m_global_int( NULL )
        , m_global_un( NULL )
        , m_tech_dup(0) {
    }

    pairwise_statistic( sequence_type & b, accum_type & d, accum_type & i, accum_type & u, double w = 1.0 ) :
        _base( &b )
        , base_weight(w)
        , m_global_diff( &d )
        , m_global_int( &i )
        , m_global_un( &u )
        , m_tech_dup(0) {
    }

    void update( sequence_type & seq, double w = 1.0 ) {
        typename sequence_type::data_citerator sit = _base->begin(), sit2 = seq.begin();

        unsigned int _diff = 0, _int = 0, _un = 0;

        while( true ) {
            if( sit == _base->end() ) {
                while( sit2 != seq.end() ) {
                    ++sit2;
                    ++_diff;
                    ++_un;
                }
                break;
            }

            if( sit2 == seq.end() ) {
                while( sit != _base->end() ) {
                    ++sit;
                    ++_diff;
                    ++_un;
                }
                break;
            }

            if( *sit < *sit2 ) {
                ++sit;
                ++_diff;
                ++_un;
            } else if( *sit2 < *sit ) {
                ++sit2;
                ++_diff;
                ++_un;
            } else {
                ++_int;
                ++_un;
                ++sit;
                ++sit2;
            }
        }

        double _w = (base_weight * w);

        m_acc_diff( (double)_diff, accum::weight = _w );
        m_acc_int( (double)_int, accum::weight = _w );
        m_acc_un( (double)_un, accum::weight = _w );

        if( !_diff ) ++m_tech_dup;

        if( m_global_int ) {
            (*m_global_int)( (double) _int, accum::weight = w );
            (*m_global_diff)( (double) _diff, accum::weight = w);
            (*m_global_un)( (double) _un, accum::weight = w );
        }
    }

    double count() {
        return accum::count( m_acc_diff );
    }

    double difference_sum() {
        return accum::weighted_sum( m_acc_diff );
    }

    double difference_mean() {
        return accum::weighted_mean( m_acc_diff );
    }

    double difference_variance() {
        return accum::weighted_variance( m_acc_diff );
    }

    double intersection_sum() {
        return accum::weighted_sum( m_acc_int );
    }

    double intersection_mean() {
        return accum::weighted_mean( m_acc_int );
    }

    double intersection_variance() {
        return accum::weighted_variance( m_acc_int );
    }

    double union_sum() {
        return accum::weighted_sum( m_acc_un );
    }

    double union_mean() {
        return accum::weighted_mean( m_acc_un );
    }

    double union_variance() {
        return accum::weighted_variance( m_acc_un );
    }

    sequence_type * _base;
    double          base_weight;

    accum_type   m_acc_diff, m_acc_int, m_acc_un;

    accum_type * m_global_diff, * m_global_int, * m_global_un;
    unsigned int m_tech_dup;
};

#endif  // PAIRWISE_STATISTIC_VECTOR_HPP_
