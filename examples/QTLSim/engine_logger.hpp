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
#ifndef ENGINE_CONFIGURATION_LOGGER_HPP_
#define ENGINE_CONFIGURATION_LOGGER_HPP_

#include <cstring>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/utility/clotho_strings.hpp"

#include <iostream>

void write_engine_config( const std::string & out_path );

void write_engine_config( boost::property_tree::ptree & elog );

#endif  // ENGINE_CONFIGURATION_LOGGER_HPP_
