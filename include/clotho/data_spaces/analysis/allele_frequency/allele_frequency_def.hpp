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
#ifndef CLOTHO_ALLELE_FREQUENCY_DEF_HPP_
#define CLOTHO_ALLELE_FREQUENCY_DEF_HPP_

#include <boost/property_tree/ptree.hpp>

namespace clotho {
namespace genetics {

template < class SequenceSpace, class SizeType >
class allele_frequency {
public:
    typedef SequenceSpace           space_type;
    typedef SizeType                size_type;
    typedef std::vector< size_type >   result_type;

    void evaluate( space_type & ss );
    result_type & getResults();
    void recordResults( boost::property_tree::ptree & log );

    virtual ~allele_frequency() {}
};

}   // namespace genetics
}   // namespace clotho

#endif  // CLOTHO_ALLELE_FREQUENCY_DEF_HPP_
