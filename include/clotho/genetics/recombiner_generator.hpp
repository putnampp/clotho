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
#ifndef RECOMBINER_GENERATOR_HPP_
#define RECOMBINER_GENERATOR_HPP_

#include "clotho/utility/random_generator.hpp"

//#include "recombiner.hpp"

namespace clotho {
namespace utility {

template < class URNG, class Sequence, class Engine >
class random_generator< URNG, recombiner< Sequence, Engine > > {
public:
    typedef random_generator< URNG, recombiner< Sequence, Engine > >    self_type;
    typedef URNG                                                        rng_type;
    typedef recombiner< Sequence, Engine >                              result_type;

    typedef clotho::utility::random_generator< URNG, Engine >           engine_generator_type;
    typedef typename engine_generator_type::result_type                 engine_type;

    random_generator( URNG & rng, boost::property_tree::ptree & config ) :
        m_rng( &rng )
        , m_eng_gen( rng, config ) {
    }

    random_generator( URNG & rng, engine_generator_type & eng_gen ) :
        m_rng( &rng )
        , m_eng_gen( eng_gen ) {
    }

    result_type operator()( ) {
        engine_type eng = m_eng_gen();
        return result_type( eng );
    }

protected:
    rng_type                * m_rng;
    engine_generator_type   m_eng_gen;
};

}   // namespace utility {
}   // namespace clotho {

#endif  // RECOMBINER_GENERATOR_HPP_
