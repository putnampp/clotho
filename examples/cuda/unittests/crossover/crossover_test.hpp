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
#ifndef CROSSOVER_TEST_HPP_
#define CROSSOVER_TEST_HPP_


template < class CrossType >
struct crossover_test {
    typedef CrossType crossover_type;

    typedef thrust::device_vector< typename CrossType::real_type >          random_vector;
    typedef thrust::device_vector< typename CrossType::allele_type >        allele_vector;
    typedef thrust::device_vector< typename CrossType::event_count_type >   count_vector;
    typedef thrust::device_vector< typename CrossType::int_type >           sequence_vector;
    typedef std::vector< typename CrossType::real_type >                    host_vector; 

    typedef typename random_vector::iterator    random_iterator;
    typedef typename count_vector::iterator     count_iterator;
    typedef typename allele_vector::iterator    allele_iterator;
    typedef typename sequence_vector::iterator  sequence_iterator;
    typedef typename host_vector::iterator      host_iterator;

    random_vector       rand_pool;
    allele_vector       allele_list;
    count_vector        event_list;
    sequence_vector     sequences;

    crossover_type      ct;

    template < class AlleleGenerator >
    void initialize( AlleleGenerator & aGen, size_t N ) {
        unsigned int tail = N % crossover_type::comp_cap_type::THREADS_PER_BLOCK;
        if( tail ) {
            N += (crossover_type::comp_cap_type::THREADS_PER_BLOCK - tail);
            std::cerr << "Warning: Increased Allele Size to - " << N << std::endl;
        }

        allele_list.resize(N);
        aGen( allele_list, N );
    }

    template < class CountGenerator, class EventGenerator >
    void simulate( CountGenerator & cGen, EventGenerator & eGen, size_t N ) {
        unsigned int nEvents = cGen( event_list, N );
        std::cout << "Events generated: " << nEvents << std::endl;

        rand_pool.resize( nEvents + N );

        eGen( rand_pool );

        unsigned int sequence_width = allele_list.size() / crossover_type::ALLELE_PER_INT;

        sequences.resize( N * sequence_width );

        typename crossover_type::real_type * pool = rand_pool.data().get();
        typename crossover_type::allele_type * alleles = allele_list.data().get();
        typename crossover_type::event_count_type * events = event_list.data().get();
        typename crossover_type::int_type * seqs = sequences.data().get();

        ct( pool, alleles, events, seqs, N, allele_list.size(), sequence_width );

        cudaDeviceSynchronize();
    }

    void get_state( boost::property_tree::ptree & s ) {
        ct.get_state( s );
    }
};

#endif  // CROSSOVER_TEST_HPP_
