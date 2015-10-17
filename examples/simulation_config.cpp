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
#include "simulation_config.h"
#include "common_commandline.h"

#include <sstream>
#include <boost/algorithm/string.hpp>

#include "generation_parameter.hpp"
#include "population_parameter.hpp"
#include "logging_parameter.hpp"
#include "clotho/random/seed_parameter.hpp"
#include "clotho/mutation/mutation_rate_parameter.hpp"
#include "clotho/recombination/recombination_rate_parameter.hpp"

#include "simulation_config_json.h" // json keys

#define RESET_SS( x ) x.str(""); x.clear();

/*
void parse_config( simulation_config & cfg ) {
    assert( boost::algorithm::iends_with( cfg.cfg_path, ".json") );

    boost::property_tree::ptree jfile;

    boost::property_tree::read_json(cfg.cfg_path, jfile);

    assert( jfile.get_child_optional( CONFIG_BLOCK_K ) != boost::none );

    std::ostringstream oss;

    oss << CONFIG_BLOCK_K << "." << GEN_BLOCK_K << "." << SIZE_K;
    cfg.nGen = jfile.get<unsigned int>( oss.str(), DEFAULT_GENERATIONS );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << POP_BLOCK_K << "." << SIZE_K;
    cfg.nPop = jfile.get< unsigned int >( oss.str(), DEFAULT_POPULATION_SIZE );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << RNG_BLOCK_K << "." << SEED_K;

    cfg.seed = jfile.get< unsigned int>(oss.str(), DEFAULT_SEED);

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << REC_BLOCK_K << "." << RATE_PER_REGION_K;
    cfg.rho = jfile.get< double >( oss.str(), DEFAULT_RECOMB_RATE );

    RESET_SS(oss)

    oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
    cfg.mu = jfile.get< double >( oss.str(), DEFAULT_MUTATION_RATE );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << REPETITION_K;
    cfg.nRep = jfile.get< unsigned int >( oss.str(), 1 );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << LOG_BLOCK_K << "." << PERIOD_K;
    cfg.log_period = jfile.get< unsigned int >( oss.str(), -1);
}*/

void parse_config( simulation_config & cfg ) {
    assert( boost::algorithm::iends_with( cfg.cfg_path, ".json") );

    typedef double real_type;

    boost::property_tree::ptree jfile;
    boost::property_tree::read_json(cfg.cfg_path, jfile);

    boost::property_tree::ptree config;
    config = jfile.get_child( CONFIG_BLOCK_K, config );

    cfg.nRep = config.get< unsigned int >( REPETITION_K, 1 );

    generation_parameter gen_param( config );
    cfg.nGen = gen_param.m_size;

    population_parameter pop_param( config );
    cfg.nPop = pop_param.m_size;

    mutation_rate_parameter< real_type > mut_param( config );
    cfg.mu = mut_param.m_mu;

    recombination_rate_parameter< real_type > rec_param( config );
    cfg.rho = rec_param.m_rho;

    logging_parameter log_param( config );
    cfg.log_period = log_param.m_period;

    seed_parameter<typename simulation_config::seed_type > seed_param( config );
    cfg.seed = seed_param.m_seed;
}


void add_config( boost::property_tree::ptree & log, const simulation_config & sim ) {
    typedef double real_type;

    std::ostringstream oss;

    // NOTE: config_file and output_file are expected to be command line parameters
//    oss  << CONFIG_K;
//    log.put( oss.str(), sim.cfg_path );
    log.put( CONFIG_K, sim.cfg_path );

//    RESET_SS( oss )
//    oss << OUTPUT_K;
//    log.put( oss.str(), sim.out_path );
    log.put( PREFIX_K, sim.out_path );

//    RESET_SS(oss)
//    oss << REPETITION_K;
//    log.put( oss.str(), sim.nRep );
    log.put( REPETITION_K, sim.nRep );

//    RESET_SS(oss)
//    oss << LOG_BLOCK_K << "." << PERIOD_K;
//    log.put( oss.str(), sim.log_period);
    logging_parameter log_param( sim.log_period );
    log_param.write_parameter( log );

//    RESET_SS( oss )
//    oss << RNG_BLOCK_K << "." << SEED_K;
//    log.put( oss.str(), sim.seed );
    seed_parameter< typename simulation_config::seed_type > seed_param( sim.seed );
    seed_param.write_parameter( log );

//    RESET_SS( oss )
//    oss << GEN_BLOCK_K << "." << SIZE_K;
//    log.put( oss.str(), sim.nGen );
    generation_parameter gen_param( sim.nGen );
    gen_param.write_parameter( log );

//    RESET_SS( oss )
//    oss << POP_BLOCK_K << "." << SIZE_K;
//    log.put( oss.str(), sim.nPop );
    population_parameter pop_param( sim.nPop );
    pop_param.write_parameter( log );

    mutation_rate_parameter< real_type > mut_param( sim.mu );
    mut_param.write_parameter( log );

    recombination_rate_parameter< real_type > rec_param( sim.rho );
    rec_param.write_parameter( log );

//    RESET_SS( oss )
//    oss << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
//    log.put( oss.str(), sim.mu );
//
//    RESET_SS( oss )
//    oss << MUT_BLOCK_K << "." << RATE_PER_BASE_K;
//    log.put( oss.str(), 0.00000001 );
//
//    RESET_SS( oss )
//    oss << REC_BLOCK_K << "." << RATE_PER_REGION_K;
//    log.put( oss.str(), sim.rho );
//
//    RESET_SS( oss )
//    oss << REC_BLOCK_K << "." << RATE_PER_BASE_K;
//    log.put( oss.str(), 0.00000001 );

    RESET_SS(oss)
    oss << OPT_BLOCK_K << "." << CHECK_SELECTED_K;
#ifdef CHECK_SELECTED
    log.put( oss.str(), true );
#else
    log.put( oss.str(), false );
#endif
}
