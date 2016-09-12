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
#ifndef QTLSIM_CONFIG_HPP_
#define QTLSIM_CONFIG_HPP_

#include "../config.h"

#include "qtlsim_logger.hpp"
#include "log_writer.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "clotho/configuration_manager/config_manager.hpp"

#include "../generation_parameter.hpp"
#include "../population_parameter.hpp"
#include "options/configuration_option.hpp"
#include "options/log_prefix_option.hpp"

#include "clotho/genetics/population_growth_toolkit.hpp"

#include "clotho/utility/timer.hpp"
#include "clotho/utility/log_helper.hpp"
#include "clotho/utility/clotho_strings.hpp"

#include "qtlsim_common.hpp"

// typedefs
typedef clotho::configuration_manager::config_manager       config_manager_type;

#ifdef USE_DOUBLE_REAL
typedef double          real_type;
#else
typedef float           real_type;
#endif
typedef unsigned int    int_type;

typedef clotho::utility::timer                              timer_type;

#ifdef USE_ROW_ALIGNMENT
#define __ALIGNMENT_TYPE__ clotho::genetics::row_grouped< 1 >
#elif USE_ROW_VECTOR
#define __ALIGNMENT_TYPE__ clotho::genetics::row_vector
#else
#define __ALIGNMENT_TYPE__ clotho::genetics::column_aligned
#endif  // USE_ROW_ALIGNMENT

#endif  // QTLSIM_CONFIG_HPP_
