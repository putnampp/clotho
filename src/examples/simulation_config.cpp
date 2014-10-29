#include "simulation_config.h"

#include <sstream>
#include <boost/algorithm/string.hpp>

#define RESET_SS( x ) x.str(""); x.clear();

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
    cfg.mu = jfile.get< double >( oss.str(), DEFAULT_RECOMB_RATE );

    RESET_SS(oss)

    oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
    cfg.mu = jfile.get< double >( oss.str(), DEFAULT_MUTATION_RATE );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << REPETITION_K;
    cfg.nRep = jfile.get< unsigned int >( oss.str(), 1 );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << LOG_FREQUENCY_K;
    cfg.log_period = jfile.get< unsigned int >( oss.str(), -1);
}

void add_config( boost::property_tree::ptree & log, const simulation_config & sim ) {
    std::ostringstream oss;

    // NOTE: config_file and output_file are expected to be command line parameters
    oss << CONFIG_BLOCK_K << ".config_file";
    log.put( oss.str(), sim.cfg_path );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << ".output_path";
    log.put( oss.str(), sim.out_path );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << REPETITION_K;
    log.put( oss.str(), sim.nRep );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << LOG_FREQUENCY_K;
    log.put( oss.str(), sim.log_period);

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << RNG_BLOCK_K << "." << SEED_K;
    log.put( oss.str(), sim.seed );
    
    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << GEN_BLOCK_K << "." << SIZE_K;
    log.put( oss.str(), sim.nGen );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << POP_BLOCK_K << "." << SIZE_K;
    log.put( oss.str(), sim.nPop );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_REGION_K;
    log.put( oss.str(), sim.mu );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << MUT_BLOCK_K << "." << RATE_PER_BASE_K;
    log.put( oss.str(), 0.00000001 );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << REC_BLOCK_K << "." << RATE_PER_REGION_K;
    log.put( oss.str(), sim.rho );

    RESET_SS( oss )
    oss << CONFIG_BLOCK_K << "." << REC_BLOCK_K << "." << RATE_PER_BASE_K;
    log.put( oss.str(), 0.00000001 );

    RESET_SS(oss)
    oss << CONFIG_BLOCK_K << "." << OPT_BLOCK_K << "." << CHECK_SELECTED_K;
#ifdef CHECK_SELECTED
    log.put( oss.str(), true );
#else
    log.put( oss.str(), false );
#endif
}
