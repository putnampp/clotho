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
#include "population_graphs.h"

#include <sstream>
#include <string>

void logMutationFrequencies( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, std::map< size_t, size_t > & freq_dist, typename locus_bitset::alphabet_t::pointer alpha ) {
    boost::property_tree::ptree v, _v0, _v1;
    _v0.put("", "Allele");
    _v1.put("", "Frequency");
    v.push_back( std::make_pair("", _v0));
    v.push_back( std::make_pair("", _v1));
    p.add_child("frequency.y.smps", v);

    boost::property_tree::ptree x, _x0, _x1;
    _x0.put("", "Region relative position of Allele" );
    _x1.put("", "Frequency of allele in population" );
    x.push_back( std::make_pair("", _x0 ) );
    x.push_back( std::make_pair("", _x1 ) );
    p.add_child( "frequency.x.Description", x);

    boost::property_tree::ptree d,s;

    locus_bitset::alphabet_t::active_iterator alpha_it = alpha->active_begin();
    for( std::vector< size_t >::iterator it = frequencies.begin(); it != frequencies.end(); ++it ) {
        if( *it > 0 ) {
            assert( alpha->checkFreeStatus( it - frequencies.begin()) );
            boost::property_tree::ptree x,y,z, _s;
            x.put( "", alpha_it->first);
            y.put( "", (*it));

            z.push_back( std::make_pair("", x ));
            z.push_back( std::make_pair("", y ));

            std::ostringstream oss;
            oss << (it - frequencies.begin() );

            _s.put("", oss.str());
            s.push_back( std::make_pair("", _s));
            d.push_back( std::make_pair("", z ));

            std::map< size_t, size_t >::iterator it2 = freq_dist.find((*it));
            if( it2 == freq_dist.end() ) {
                freq_dist.insert( std::make_pair( (*it), 1 ) );
            } else {
                ++it2->second;
            }
        }
        ++alpha_it;
    }

    p.add_child( "frequency.y.vars", s );
    p.add_child( "frequency.y.data", d );

    boost::property_tree::ptree graph_opts;
    graph_opts.put("graphType", "Scatter2D");
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Allele");
        tmp.push_back( std::make_pair("", t ) );
        graph_opts.add_child("xAxis", tmp);
    }
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Frequency");
        tmp.push_back( std::make_pair( "",t));
        graph_opts.add_child("yAxis", tmp);
    }
    graph_opts.put( "title", "Allele Frequency" );
    p.add_child( "frequency.graph_opts", graph_opts );
}

void logMutationBinFrequencies( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, unsigned int bin_count, typename locus_bitset::alphabet_t::pointer alpha ) {
    boost::property_tree::ptree v, _v0, _v1, _v2;
    std::string key = "genome_bin";
    _v0.put("", "Bin");
    _v1.put("", "Offset");
    _v2.put("", "Frequency");
    v.push_back( std::make_pair("", _v0));
    v.push_back( std::make_pair("", _v1));
    v.push_back( std::make_pair("", _v2));
    p.add_child( key +".y.smps", v);

    boost::property_tree::ptree x, _x0, _x1, _x2;
    _x0.put("", "Bin index" );
    _x1.put("", "Genomic offset relative to bin" );
    _x2.put("", "Frequency of allele in population" );
    x.push_back( std::make_pair("", _x0 ) );
    x.push_back( std::make_pair("", _x1 ) );
    x.push_back( std::make_pair("", _x2 ) );

    p.add_child( key + ".x.Description", x);

    boost::property_tree::ptree d,s;

    typename locus_bitset::alphabet_t::active_iterator alpha_it = alpha->active_begin();

    typedef  std::vector< std::map< double, size_t > > bin_freq_type;
    bin_freq_type bin_freq( bin_count, std::map<double, size_t>() );

    for( std::vector< size_t >::iterator it = frequencies.begin(); it != frequencies.end(); ++it ) {
        if( *it > 0 ) {
            assert( alpha->checkFreeStatus( it - frequencies.begin()) );
            unsigned int bin_idx = alpha_it->first * bin_count;

            double lo = (double)(bin_idx) / (double) bin_count;

            double offset = alpha_it->first - lo;

            bin_freq[bin_idx].insert( std::make_pair( offset, (*it)));
        }
        ++alpha_it;
    }

    unsigned int bin_idx = 0;
    for( bin_freq_type::iterator it = bin_freq.begin(); it != bin_freq.end(); ++it ) {
        double lo = (double)(bin_idx) / (double) bin_count;
        boost::property_tree::ptree w;
        w.put( "", lo);
        unsigned int offset = 0;
        for( std::map< double, size_t >::iterator it2 = it->begin(); it2 != it->end(); ++it2 ) {
            boost::property_tree::ptree x,y,z, _s;
            x.put( "", it2->first);
            y.put( "", it2->second);
            z.push_back( std::make_pair("", w ));
            z.push_back( std::make_pair("", x ));
            z.push_back( std::make_pair("", y ));

            d.push_back( std::make_pair("", z ));

            std::ostringstream oss;
            oss << bin_idx << "_" << offset++;
            _s.put("", oss.str());
            s.push_back( std::make_pair("", _s));
        }
        ++bin_idx;
    }

    p.add_child( key + ".y.vars", s );
    p.add_child( key + ".y.data", d );

    boost::property_tree::ptree graph_opts;

    graph_opts.put("graphType", "Scatter3D");

    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Bin");
        tmp.push_back( std::make_pair("", t ) );
        graph_opts.add_child("xAxis", tmp);
    }
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Offset");
        tmp.push_back( std::make_pair( "",t));
        graph_opts.add_child("zAxis", tmp);
    }
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Frequency");
        tmp.push_back( std::make_pair("", t ));
        graph_opts.add_child("yAxis", tmp);
    }
    graph_opts.put( "title", "Mutation Distribution grouped by bin" );

    p.add_child( key + ".graph_opts", graph_opts );
}

void logMutationDistribution( boost::property_tree::ptree & p, std::map< size_t, size_t > & freq_dist ) {
    boost::property_tree::ptree v, _v0, _v1;
    _v0.put("", "Frequency");
    _v1.put("", "Count");
    v.push_back( std::make_pair("", _v0));
    v.push_back( std::make_pair("", _v1));
    p.add_child("distribution.y.smps", v);

    boost::property_tree::ptree x, _x0, _x1;
    _x0.put("", "Region relative position of Allele" );
    _x1.put("", "Frequency of allele in population" );
    x.push_back( std::make_pair("", _x0 ) );
    x.push_back( std::make_pair("", _x1 ) );
    p.add_child( "distribution.x.Description", x);
    boost::property_tree::ptree d,s;

    unsigned int i = 0;
    for( std::map< size_t, size_t >::iterator it = freq_dist.begin(); it != freq_dist.end(); ++it ) {
        boost::property_tree::ptree x,y,z, _s;
        x.put( "", it->first);
        y.put( "", it->second);

        z.push_back( std::make_pair("", x ) );
        z.push_back( std::make_pair("", y ) );

        std::ostringstream oss;
        oss << i++;

        _s.put("", oss.str());
        s.push_back( std::make_pair("", _s));
        d.push_back( std::make_pair("", z ));
    }

    p.add_child( "distribution.y.vars", s );
    p.add_child( "distribution.y.data", d );

    boost::property_tree::ptree graph_opts;
    graph_opts.put("graphType", "Scatter2D");
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Frequency");
        tmp.push_back( std::make_pair("", t ) );
        graph_opts.add_child("xAxis", tmp);
    }
    {
        boost::property_tree::ptree tmp, t;
        t.put("", "Count");
        tmp.push_back( std::make_pair( "",t));
        graph_opts.add_child("yAxis", tmp);
    }
    graph_opts.put( "title", "Allele Frequency Distribution" );
    p.add_child( "distribution.graph_opts", graph_opts );
}

void logMutationStats( boost::property_tree::ptree & p, std::vector< size_t > & frequencies, typename locus_bitset::alphabet_t::pointer alpha ) {
    std::map< size_t, size_t > freq_dist;
    logMutationFrequencies( p, frequencies, freq_dist, alpha );
    logMutationBinFrequencies( p, frequencies, 1000, alpha );
    logMutationDistribution( p, freq_dist );
}
