import json, itertools, os, random

def make_generation( g ):
    return { "size": "{0}".format(g) }

def make_log_sample( s, pw ):
    if pw:
        pw = "true"
    else:
        pw = "false"

    return { "size": "{0}".format(s), "pairwise": "{0}".format(pw) }

def make_log( p, samples=None ):
    ret = {"period": "{0}".format(p) }
    if samples != None:
        ret[ "sampling_size" ] = [ make_log_sample( s, pw ) for s, pw in samples ]

    return ret

def make_random( seed=None ):
    if seed == None:
        return { "seed": "" }
    else:
        return { "seed": "{0}".format(seed) }

def make_multithread( t ):
    return { "T": "{0}".format(t) }

def make_normal_traits( N, mean="0", sigma="1" ):
    return { "N": "{0}".format(N), "mean": "{0}".format(0), "sigma": "{0}".format(sigma) }

def make_allele_neutrality( p ):
    return { "p": "{0}".format(p) }

def make_recombination( rho, bias=0.5 ):
    return { "sequence_bias": { "p": "{0}".format(bias) }, "rho": "{0}".format(rho) }

def make_mutation( mu ):
    return { "mu": "{0}".format(mu) }

def make_pop_growth( method, params ):
    return { "name": "{0}".format(method), "params": params }

def make_constant_population( s ):
    return make_pop_growth( "linear", { "A": "0", "B": "{0}".format(s) } )

def make_fitness( method, params ):
    return { "name": "{0}".format(method), "params": params }

def make_quadratic_fitness( scale ):
    return make_fitness( "quadratic", { "scale": "{0}".format(scale) } )

def make_initial_allele_frequency_distribution( A ):
    ad = []
    while len(ad) < A:
        ad.append( random.uniform(0.001, 1.0) )

    return { "allele_distribution" : ad }


def make_constant_quadratic_configuration( G, P, neut, mu, rho, thread, A, quad_fit_scale ):
    pconfig = {}
    pconfig[ "generations" ] = make_generation( G )
    pconfig[ "population_growth" ] = make_constant_population( P )
    pconfig[ "neutral" ] = make_allele_neutrality( neut )
    pconfig[ "mutation" ] = make_mutation( mu )
    pconfig[ "recombination" ] = make_recombination( rho )
    pconfig[ "fitness" ] = make_quadratic_fitness( quad_fit_scale )
    pconfig[ "traits" ] = make_normal_traits( "1", "0.0", "1.0" )
    pconfig[ "random_number" ] = make_random()
    pconfig[ "multithreading" ] = make_multithread( thread )
    pconfig[ "log" ] = make_log( 100 )

    if A > 0:
        pconfig[ "initialize" ] = make_initial_allele_frequency_distribution( A )

    return pconfig
    

def recombination_performance_experiment( path_prefix ):
    # generations
    G = [ "100" ] 
    # population size
    P = [ "10000" ]

    p_neutrality = [ "0.0", "1.0" ]
    mu = [ "0.00000001" ]
    rho = [ "1.0", "10.0", "23.0", "35.0", "50.0", "100.0" ]
    threads = [ "0", "11" ]
    init_alleles = [ 500, 5000, 50000, 100000 ]
    # quadratic scale
    scale = [ "4.0" ]

    if path_prefix != "":
        path_prefix = "{0}/".format(path_prefix)

    for g, p, n, m, r, t, a, s in itertools.product(G, P, p_neutrality, mu, rho, threads, init_alleles, scale ):
        print("{4};{5};{0} x {1} x {2} x {3}".format( n, m, r, t, g, p))
        config = make_constant_quadratic_configuration( g, p, n, m, r, t, a, s )

        path="{0}neutral={1}/mu={2}/rho={3}/A={4}/T={5}".format( path_prefix, n, m, r, a, t )
        if not os.path.exists( path ):
            os.makedirs( path )

        with open("{0}/config.json".format(path), "w") as of:
            json.dump( config, of, indent=2 )

p=""
recombination_performance_experiment(p)
