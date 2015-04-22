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
#ifndef IGENERATOR_HPP_
#define IGENERATOR_HPP_

template < class Element >
struct igenerator {
    typedef Element result_type;

    virtual const std::string name() const = 0;
    virtual result_type operator()() = 0;
};

template < class URNG, class Element >
struct iran_generator : virtual public igenerator< Element > {

    virtual URNG *  getRNG() = 0;
};

#endif  // IGENERATOR_HPP_
