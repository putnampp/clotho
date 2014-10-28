#ifndef SIMULATION_CONFIG_H_
#define SIMULATION_CONFIG_H_

#include "config.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "simulation_config_json.h"

struct simulation_config {
    unsigned int    nGen, nPop;
    double          mu, rho;

    unsigned int    seed;

    string          cfg_path, out_path;
    simulation_config( unsigned int gen = 0, unsigned int pop = 0, double m = 0.0, double r = 0.0, unsigned int s = 0, const string & cpath = "", const string & opath = "" ) :
        nGen(gen)
        , nPop(pop)
        , mu(m)
        , rho(r)
        , seed(s)
        , cfg_path(cpath)
        , out_path(opath)
    {}

    simulation_config( const simulation_config & s) :
        nGen(s.nGen)
        , nPop(s.nPop)
        , mu(s.mu)
        , rho(s.rho)
        , seed(s.seed)
        , cfg_path(s.cfg_path)
        , out_path(s.out_path)
    {}
};

void parse_config( simulation_config & cfg );
void add_config( boost::property_tree::ptree & log, const simulation_config & sim );

#endif  // SIMULATION_CONFIG_H_
