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
#ifndef SIMULATION_CONFIG_H_
#define SIMULATION_CONFIG_H_

#include "config.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//#include "simulation_config_json.h"

struct simulation_config {
    typedef unsigned int seed_type;

    unsigned int    nGen, nPop, nRep;
    double          mu, rho;

    seed_type    seed;

    string          cfg_path, out_path;
    unsigned int    log_period;
    simulation_config( unsigned int gen = 0, unsigned int pop = 0, unsigned int rep = 0, double m = 0.0, double r = 0.0, unsigned int s = 0, const string & cpath = "", const string & opath = "", unsigned int log = 0 ) :
        nGen(gen)
        , nPop(pop)
        , nRep(rep)
        , mu(m)
        , rho(r)
        , seed(s)
        , cfg_path(cpath)
        , out_path(opath)
        , log_period(log) {
    }

    simulation_config( const simulation_config & s) :
        nGen(s.nGen)
        , nPop(s.nPop)
        , nRep(s.nRep)
        , mu(s.mu)
        , rho(s.rho)
        , seed(s.seed)
        , cfg_path(s.cfg_path)
        , out_path(s.out_path)
        , log_period(s.log_period) {
    }
};

void parse_config( simulation_config & cfg );
void add_config( boost::property_tree::ptree & log, const simulation_config & sim );

#endif  // SIMULATION_CONFIG_H_
