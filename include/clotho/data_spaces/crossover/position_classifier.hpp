//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef CLOTHO_POSITION_CLASSIFIER_HPP_
#define CLOTHO_POSITION_CLASSIFIER_HPP_

#include "clotho/classifiers/binary_classifier.hpp"
#include "clotho/classifiers/region_classifier.hpp"
#include "clotho/classifiers/tags/is_odd_tag.hpp"

namespace clotho {
namespace genetics {

template < class PositionType >
class PositionClassifier {
public:
    typedef PositionType    position_type;

    typedef clotho::classifiers::region_classifier< PositionType > base_classifier_type;
    typedef clotho::classifiers::binary_classifier< base_classifier_type, clotho::classifiers::tags::is_odd_tag > classifier_type;

    typedef typename base_classifier_type::param_type event_type;

    typedef PositionClassifier< PositionType >  self_type;

    PositionClassifier( position_type * pos, const event_type & events ) :
        m_positions( pos )
        , m_cfier( events )
    {}

    PositionClassifier( const self_type & other ) : 
        m_positions( other.m_positions )
        , m_cfier( other.m_cfier )
    {}

    bool operator()( unsigned int offset ) {
//#ifdef DEBUGGING
//        BOOST_LOG_TRIVIAL( debug ) << "Comparing offset: " << offset << " (" << m_positions[offset] << ")";
//#endif  // DEBUGGING

        return m_cfier( m_positions[ offset ] );
    }

    size_t event_count() const {
        return m_cfier.event_count();
    }

    virtual ~PositionClassifier() {}
protected:
    position_type   * m_positions;
    classifier_type m_cfier;
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_POSITION_CLASSIFIER_HPP_
