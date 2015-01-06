#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include <boost/algorithm/string.hpp>

using std::string;
using std::ostream;
using std::ofstream;

namespace po=boost::program_options;

const string HELP_K = "help";
const string CONFIG_K = "config";
const string INPUT_K = "input";
const string OUTPUT_K = "output";

const string OPERATION_K = "operation";
const string SELECT_K = "select";
const string SPLIT_K = "split";
const string TRANSPOSE_K = "transpose";
const string COMBINE_K = "combine";

const string FORMAT_K = "format";
const string PREFIX_K = "prefix";
const string VALUE_K = "value";
const string KEYS_K = "keys";

const string JSON_K = "json";
const string GNUPLOT_K = "gnuplot";

const string OUTPUT_FORMAT_K = OUTPUT_K + "." + FORMAT_K;
const string OUTPUT_PREFIX_K = OUTPUT_K + "." + PREFIX_K;

const string INPUT_VALUE_K = INPUT_K + "." + VALUE_K;

int parse_configuration( int argc, char ** argv, boost::property_tree::ptree & cfg );
int command_line( int argc, char ** argv, po::variables_map & vm );
int transform( po::variables_map & vm, boost::property_tree::ptree & p );

void parse_data( boost::property_tree::ptree & cfg );

typedef void (*parser_op)( boost::property_tree::ptree & parameters );

// select keys from each input file into a
void select_op( boost::property_tree::ptree & parameters );
void select_op_gnuplot( boost::property_tree::ptree & params );

// combine all keys from each input file into a single output file
void combine_op( boost::property_tree::ptree & parameters );
void combine_op_json( boost::property_tree::ptree & parameters );
void transpose_op( boost::property_tree::ptree & params );

void write_json( const std::string & path, boost::property_tree::ptree & results );
void write_gnuplot( const std::string & path, const std::string & header, boost::property_tree::ptree & results );

int main( int argc, char ** argv ) {

    boost::property_tree::ptree cfg;
    int res = 0;
    if( (res = parse_configuration( argc, argv, cfg )) ) {
        return res;
    }

    BOOST_FOREACH( auto& c, cfg ) {
        std::cerr << c.first << std::endl;
        parse_data( c.second );
    }
}

void parse_data( boost::property_tree::ptree & cfg ) {
    string op = cfg.get< string >( OPERATION_K, "");

    if( op.empty() ) return;

    if( op == SELECT_K ) {
        select_op( cfg );
    } else if ( op == COMBINE_K ) {
//        combine_op( cfg );
    } else if( op == TRANSPOSE_K ) {
//        transpose_op( cfg );
    }
}

int parse_configuration( int argc, char ** argv, boost::property_tree::ptree & cfg ) {
    po::variables_map vm;
    int res = command_line( argc, argv, vm );

    if( res ) {
        return res;
    }

    string config_file = vm[ CONFIG_K ].as< string >();
    if( config_file.empty() ) {
        res = transform( vm, cfg );
    } else {
        boost::property_tree::read_json( config_file, cfg );
    }

    return res;
}


int command_line( int argc, char ** argv, po::variables_map & vm ) {
    po::options_description gen( "General" );
    gen.add_options()
    ((HELP_K + ",h").c_str(), "Print this" )
    ;

    po::options_description io_params("I/O Parameters");
    io_params.add_options()
    ( (CONFIG_K + ",c" ).c_str(), po::value< string >()->default_value(""), "Configuration file containing files to be parsed" )
    ;

    po::options_description cmdline;
    cmdline.add( gen ).add( io_params );
    po::store( po::command_line_parser( argc, argv ).options(cmdline).run(), vm );

    if( vm.count( HELP_K ) ) {
        std::cout << cmdline << std::endl;
        return 1;
    }

    return 0;
}

int transform( po::variables_map & vm, boost::property_tree::ptree & p ) {
    return 0;
}

void select_op( boost::property_tree::ptree & params ) {

    string frmt = params.get< string >( OUTPUT_FORMAT_K, "" );
    string output = params.get< string >( OUTPUT_PREFIX_K, "" );

    std::vector< string > keys;
    BOOST_FOREACH( auto& v, params.get_child( KEYS_K ) ) {
        std::ostringstream oss;
        oss << v.second.data();
        keys.push_back( oss.str() );
    }

    BOOST_FOREACH( auto& v, params.get_child( INPUT_K ) ) {
        std::ostringstream oss;
        oss << v.second.data();

        string path = oss.str();

        if( boost::algorithm::iends_with(path, "." + JSON_K) ) {
            std::cerr << "Parsing input JSON file: " << path << std::endl;

            string tmp_path = path.substr(0, path.size() - (JSON_K.size() + 1));
            boost::property_tree::ptree infile;
            boost::property_tree::read_json( path, infile );

            BOOST_FOREACH( auto& k, keys ) {
                if( infile.get_child_optional( k ) != boost::none ) {
//                    string tmp = tmp_path + "." + k;
//                    results.add_child( k, infile.get_child( k ) );
//

                    if( boost::iequals( frmt, JSON_K ) ) {
                        write_json( output, infile.get_child(k) );
                    } else if( boost::iequals( frmt, GNUPLOT_K ) ) {
                        write_gnuplot( output, k, infile.get_child(k) );
                    }
                }
            }
        }
    }
}
/*
void combine_op_json( boost::property_tree::ptree & params ) {
    boost::property_tree::ptree results;
    assert( params.get_child_optional( INPUT_K ) != boost::none );
    assert( params.get_child_optional( KEYS_K ) != boost::none );

    std::vector< string > keys;

    BOOST_FOREACH( auto& v, params.get_child( KEYS_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();
            keys.push_back( oss.str() );
    }

    BOOST_FOREACH( auto& v, params.get_child( INPUT_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();

            string path = oss.str();

            if( boost::algorithm::iends_with(path, "." + JSON_K) ) {
                std::cerr << "Parsing input JSON file: " << path << std::endl;

                string tmp_path = path.substr(0, path.size() - (JSON_K.size() + 1));
                boost::property_tree::ptree infile;
                boost::property_tree::read_json( path, infile );

                BOOST_FOREACH( auto& k, keys ) {
                    if( infile.get_child_optional( k ) != boost::none ) {
                        string tmp = tmp_path + "." + k;
                        results.add_child( tmp, infile.get_child( k ) );
                    }
                }
            }
        }

        if( output.empty() ) {
            boost::property_tree::write_json(std::cout, results );
        } else {
            boost::property_tree::write_json( output + "." + JSON_K, results );
        }

}

void select_op_gnuplot( boost::property_tree::ptree & params ) {

}

void transpose_op( boost::property_tree::ptree & params ) {
    assert( params.get_child_optional( OUTPUT_FORMAT_K ) != boost::none );

    string frmt = params.get< string >( OUTPUT_FORMAT_K, "" );
    string output = params.get< string >( OUTPUT_PREFIX_K, "" );

    if( boost::iequals( frmt, GNUPLOT_K ) ) {
        assert( params.get_child_optional( INPUT_K ) != boost::none );
        assert( params.get_child_optional( KEYS_K ) != boost::none );

        std::vector< string > keys;

        BOOST_FOREACH( auto& v, params.get_child( KEYS_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();
            keys.push_back( oss.str() );
        }

        BOOST_FOREACH( auto& v, params.get_child( INPUT_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();

            string path = oss.str();

            typedef std::vector< string > row_type;
            typedef std::map< string, std::vector< row_type > > value_map_type;
            value_map_type data;

            if( boost::algorithm::iends_with(path, "." + JSON_K) ) {
                std::cerr << "Parsing input JSON file: " << path << std::endl;

                string tmp_path = path.substr(0, path.size() - (JSON_K.size() + 1));
                boost::property_tree::ptree infile;
                boost::property_tree::read_json( path, infile );

                BOOST_FOREACH( auto& k, keys ) {
                    if( infile.get_child_optional( k ) != boost::none ) {
//                        string tmp = tmp_path + "." + k;
                        value_map_type::iterator it = data.find( k );
                        if( it == data.end() ) {
                            std::pair< value_map_type::iterator, bool > res = data.insert( std::make_pair(k, std::vector<row_type>() ));
                            assert( res.second );
                            it = res.first;
                        }

                        unsigned int max_columns = 0;
                        unsigned int j = 0;
                        BOOST_FOREACH( auto& dval, infile.get_child(k) ) {
                            unsigned int i = 0;
                            row_type r;
                            BOOST_FOREACH( auto& pval, dval.second ) {
                                std::ostringstream oss;
                                oss << pval.second.data();
                                r.push_back( oss.str() );
                                ++i;
                            }
                            if( i > max_columns ) { max_columns = i; }

                            it->second.push_back(r);
                            ++j;
                        }
                    }
                }
            }

            if( output.empty() ) {
                BOOST_FOREACH( auto& dval, data) {
                    std::cout << "# " << dval.first << std::endl;
                    unsigned int i = 0;
                    BOOST_FOREACH( auto& r, dval.second ) {
                        std::cout << i++;
                        BOOST_FOREACH( auto& c, r ) {
                            std::cout << "," << c;
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}

void combine_op( boost::property_tree::ptree & params ) {
    assert( params.get_child_optional( OUTPUT_FORMAT_K ) != boost::none );

    string frmt = params.get< string >( OUTPUT_FORMAT_K, "" );
    string output = params.get< string >( OUTPUT_PREFIX_K, "" );

    if( boost::iequals( frmt, JSON_K ) ) {
        boost::property_tree::ptree  results;

        assert( params.get_child_optional( INPUT_K ) != boost::none );
        assert( params.get_child_optional( KEYS_K ) != boost::none );

        std::vector< string > keys;

        BOOST_FOREACH( auto& v, params.get_child( KEYS_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();
            keys.push_back( oss.str() );
        }
        BOOST_FOREACH( auto& v, params.get_child( INPUT_K ) ) {
            std::ostringstream oss;
            oss << v.second.data();

            string path = oss.str();

            if( boost::algorithm::iends_with(path, "." + JSON_K) ) {
                std::cerr << "Parsing input JSON file: " << path << std::endl;

                string tmp_path = path.substr(0, path.size() - (JSON_K.size() + 1));
                boost::property_tree::ptree infile;
                boost::property_tree::read_json( path, infile );

                BOOST_FOREACH( auto& k, keys ) {
                    if( infile.get_child_optional( k ) != boost::none ) {
                        string tmp = tmp_path + "." + k;
                        results.add_child( tmp, infile.get_child( k ) );
                    }
                }
            }
        }
    }
}*/

void write_json( const string & path, boost::property_tree::ptree & results ) {
    if( path.empty() ) {
        boost::property_tree::write_json( std::cout, results );
    } else {
        boost::property_tree::write_json( path + "." + JSON_K, results );
    }
}
/*
void write( ostream * out, boost::property_tree::ptree::value_type  & node ) {
    if( node.first.empty() ) {
        (*out) << "\n";
        BOOST_FOREACH( auto & pval, node.second ) {
            std::ostringstream oss;
            oss << pval.second.data();

            (*out) << oss.str() << ",";
        }
    } else {
        (*out) << "." << node.first;
        BOOST_FOREACH( auto& c, node.second ) {
            write( out, c );
            (*out) << "\n";
        }
    }
}*/

void write_data( ostream * out, boost::property_tree::ptree & data ) {
    unsigned int i = 0;
    BOOST_FOREACH( auto& c, data ) {
        if( c.second.empty() ) {
            std::ostringstream oss;
            oss << c.second.data();
            (*out) << "," << oss.str();
        } else {
            (*out) << ++i;
            write_data( out, c.second);
            (*out) << "\n";
        }
    }
}

void write( ostream * out, string path, boost::property_tree::ptree & tr ) {
    unsigned int i = 0;
    BOOST_FOREACH( auto& c, tr ) {
        if( c.first.empty() ) {
            (*out) << ++i;
            write_data(out, c.second);
            (*out) << "\n";
        } else {
            write( out, path + "." + c.first, c.second );
        }
    }
}

void write_gnuplot( const string & path, const string & header, boost::property_tree::ptree & results ) {
    ostream * out;
    if( path.empty() ) {
        out = &std::cout;
    } else {
        out = new ofstream(path + ".dat");
    }

    (*out) << "# " << header << "\n";
    write(out, "", results );

    if( out != &std::cout ) {
        ((ofstream *)out)->close();
    }
}
