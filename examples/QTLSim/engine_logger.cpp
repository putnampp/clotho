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
#include "engine_logger.hpp"

#define XSTR( x ) #x
#define STR( x )  XSTR( x )

const std::string ENGINE_BLOCK_K = "engine";
const std::string POWERSET_BLOCK_K = "powerset";
const std::string SUBSET_BLOCK_K = "subset";
const std::string CLASSIFIER_BLOCK_K = "classifier";

void write_engine_config( const std::string & out_path ) {

    boost::property_tree::ptree sset, recomb, pset;
    recomb.put( "tag0", STR( RECOMBINE_INSPECT_METHOD ) );
    recomb.put( "tag1", STR( BIT_WALK_METHOD ) );
    sset.put_child( "recombine", recomb );
    sset.put( TYPE_K, STR( SUBSETTYPE ) );

    pset.put_child( SUBSET_BLOCK_K, sset );
//    pset.put( "block_size", BLOCK_UNIT_SIZE );
//    compile_log.put(ENGINE_BLOCK_K + "." + REC_BLOCK_K + "." + TYPE_K, STR( (RECOMBTYPE) ) );
//    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K +"." + SIZE_K, BLOCK_UNIT_SIZE );

//    compile_log.put(ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K  +"." + TYPE_K, STR( SUBSETTYPE ) );
//    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag0", STR( RECOMBINE_INSPECT_METHOD ) );
//    compile_log.put( ENGINE_BLOCK_K + "." + POWERSET_BLOCK_K + "." + SUBSET_BLOCK_K + "." + REC_BLOCK_K + ".tag1", STR( BIT_WALK_METHOD ) );

//    compile_log.put( ENGINE_BLOCK_K + ".reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
//
//    compile_log.put( ENGINE_BLOCK_K + ".individual_selector.type", STR(IND_SELECT) );

    boost::property_tree::ptree eng;
    eng.put_child( POWERSET_BLOCK_K, pset );
    eng.put( REC_BLOCK_K + ".type", STR( (RECOMBTYPE) ) );
    eng.put( "reproduction_method.type", STR(REPRODUCTION_METHOD_TAG));
    eng.put( "individual_selector.type", STR(IND_SELECT) );
    eng.put( "description", "Simulator compiled objects; READ ONLY");

    boost::property_tree::ptree compile_log;
    compile_log.put_child( ENGINE_BLOCK_K, eng );

    if( out_path.empty() ) {
        boost::property_tree::write_json( std::cout, compile_log );
    } else {
        std::ostringstream oss;
        oss << out_path << ".engine_compilation.json";

        boost::property_tree::write_json( oss.str(), compile_log );
    }
}
