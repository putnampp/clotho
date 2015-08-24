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
#ifndef VALIDATE_EVENT_SEQUENCE_HPP_
#define VALIDATE_EVENT_SEQUENCE_HPP_

#include <vector>
#include <boost/property_tree/ptree.hpp>

template < class TableType >
bool validate_event_sequence_distribution( TableType & edt, boost::property_tree::ptree & err ) {
    int total = 0;  // total event sequences
    typename TableType::iterator it = edt.begin();

    std::vector< int > dist( edt.size(), 0 );

    typename TableType::iterator most_freq = edt.begin();

    while( it != edt.end() ) {
        total += it->second;
        dist[ it->second ] += 1;
        if( it->second > most_freq->second ) {
            most_freq = it;
        }
        ++it;
    }

    double percent_not_unique = (1.0 - ((double)dist[ 1 ] / (double) total));

    bool valid = (percent_not_unique < 0.001);  // valid if percent_unique >= 99.9% == percent_not_unique < 0.001

    if( !valid ) {
        err.put( "message", "Event sequence distribution is not sufficiently unqiue" );
        err.put( "event_sequence.not_unique.percent", percent_not_unique );
        err.put( "event_sequence.most_frequent.sequence", most_freq->first );
        err.put( "event_sequence.most_frequent.frequency", most_freq->second );
        err.put( "event_sequence.total", total );
    }

    return valid;
}
#endif  // VALIDATE_EVENT_SEQUENCE_HPP_
