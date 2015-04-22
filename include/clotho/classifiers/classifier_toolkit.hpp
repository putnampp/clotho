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
#ifndef CLASSIFIER_TOOLKIT_HPP_
#define CLASSIFIER_TOOLKIT_HPP_

#include "clotho/classifier/iclassifier_generator.hpp"

#include <unordered_map>
#include <memory>
#include <boost/property_tree/ptree.hpp>

class classifier_toolkit {
public:

    typedef std::unordered_map< std::string, iclassifier_generator * >  generator_map;
    typedef generator_map::iterator                                     generator_iterator;
    typedef generator_map::const_iterator                               generator_citerator;

    static classifier_toolkit * getInstance() {
        static classifier_toolkit inst;
        return &inst;
    }

    std::shared_ptr< iclassifier_generator > get_classifier( boost::property_tree::ptree & config );

    void classifier_configurations( boost::property_tree::ptree & config );

    void register_classifier( iclassifier_generator * cfier );

    friend std::ostream & operator<<( std::ostream & out, const classifier_toolkit & ctk);

    virtual ~classifier_toolkit();
protected:
    classifier_toolkit();

    generator_map   m_classifiers;
};

std::ostream & operator<<( std::ostream & out, const classifier_toolkit & ctk );

#endif  // CLASSIFIER_TOOLKIT_HPP_
