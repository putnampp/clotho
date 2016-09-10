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
#include "qtlsim_config.hpp"
#include "../qtl/qtl_config.hpp"

#include <iostream>

#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

//#include "clotho/data_spaces/crossover/crossover_task.hpp"
#include <boost/random/mersenne_twister.hpp>

typedef boost::random::mt19937          random_engine_type;

boost::mutex mio;

//typedef clotho::genetics::crossover_task< random_engine_type, unsigned long long, double > xover_type;
//
void my_task( int val ) {
    boost::lock_guard< boost::mutex > guard( mio );
    std::cout << "Val: " << val << std::endl;
}

void run_tasks( int start, int end, int tc = 6 ) {
    boost::asio::io_service service;
    boost::thread_group threads;

    {
        std::auto_ptr< boost::asio::io_service::work > work(new boost::asio::io_service::work( service));

        // two worker threads
        for( int i = 0; i < tc; ++i )
            threads.create_thread( boost::bind( &boost::asio::io_service::run, &service) );

        for ( int i = start; i < end; ++i ) {
            service.post( boost::bind( my_task, i ) );
        }

//        boost::this_thread::sleep_for( boost::chrono::seconds(2));
//        service.stop();
    }

    threads.join_all();
}

int main( int argc, char ** argv ) {

    po::variables_map vm;
    int ret = config_manager_type::getInstance()->parse_commandline( argc, argv, vm );
    if( ret ) return ret;

    std::string path = "sample.log";
    init_logger( path, logging::trivial::info );

    run_tasks( 0, 200 );

    run_tasks( 2000, 2600 );

    return 0;
}
