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
#include "qtlsim_logger.hpp"

void init_logger( std::string & path, logging::trivial::severity_level level ) {
    logging::add_common_attributes();

    typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;
    boost::shared_ptr< text_sink > sink = boost::make_shared< text_sink >();

    sink->locked_backend()->add_stream( boost::make_shared< std::ofstream >( path ) );

   logging::core::get()->add_sink(sink);

    auto fmtThreadID = expr::attr< attributes::current_thread_id::value_type >("ThreadID");
    auto fmtSeverity = expr::attr< logging::trivial::severity_level >("Severity");
    auto fmtTS = expr::format_date_time< boost::posix_time::ptime >( "TimeStamp", "%H:%M:%S.%f" );
    
    logging::formatter logFmt = expr::format( "%4%,%1%,%2%,%3%" ) % fmtThreadID % fmtSeverity % expr::smessage % fmtTS;

    sink->set_formatter( logFmt );
    sink->set_filter( logging::trivial::severity >= level );
}
