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
#include <boost/test/unit_test.hpp>

#include "../unittest_config.h"

#include "clotho/classifiers/region_classifier.hpp"

BOOST_AUTO_TEST_SUITE( test_classifiers )

BOOST_AUTO_TEST_CASE( region_classifier_test ) {
    clotho::classifiers::region_classifier< double > d_class;
}

BOOST_AUTO_TEST_SUITE_END()
