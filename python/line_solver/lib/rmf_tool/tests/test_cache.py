"""
Test that computes the mean field approximation and refined mean field approximation for the cache replacement policy 'RANDOM'. 
This file should test if the 'dimension reduction' works.
Compare the computed value with a value already stored in a pickle file
"""
import pickle
import numpy as np
from approximately_equal import approximately_equal

import os
PWD=os.getcwd()
if PWD[-5:] == 'tests':
    CACHE_DIR = 'output_tests'
else:
    CACHE_DIR = 'tests/output_tests'

import sys
sys.path.append('../')
sys.path.append('.')
import src.rmf_tool as rmf

def cacheModel(N, M, alpha):
    """
    N = number of items, M = cache size, alpha = parameter of the Zipf distribution.
    """
    ddpp = rmf.DDPP()

    assert M<N, "Cache size should be smaller than the number of items"

    pop = np.arange(1, N+1)**(-alpha)
    pop /= np.sum(pop)
    print(pop)

    def exchange(i, j):
        """ Add i to the cache and remove j"""
        l = np.zeros(2*N)
        l[i] = -1
        l[N+i] = 1
        l[j] += 1
        l[N+j] += -1
        return l

    # We then add the transitions : 
    for i in range(N):
        for j in range(N):
            ddpp.add_transition(exchange(i, j), eval('lambda x: {}*x[{}]*x[{}]/{}'.format(pop[i], i, N+j, M)))
    initial_state = np.zeros(2*N)
    initial_state[0:M] = 1
    initial_state[N+M:2*N]=1
    ddpp.set_initial_state(initial_state)
    pi, V, W = ddpp.meanFieldExpansionSteadyState(order=1)
    return ddpp

def generate_data():
    """
    Generate all data and store them in a pickle file
    (to be used one times when the test is initialized)
    """
    data = dict([])
    for alpha in [0.5]:
        for N in [10]:
            for M in [5, 3]:
                for order in [0, 1]:
                    ddpp = cacheModel(N, M, alpha)
                    data[(N, M, alpha, order)] = ddpp.meanFieldExpansionSteadyState(order=order)
    with open('{}/cache.pickle'.format(CACHE_DIR), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

def test_cache():
    with open('{}/cache.pickle'.format(CACHE_DIR), 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
        for key in data:
            (N, M, alpha, order) = key
            print(key)
            ddpp = cacheModel(N, M, alpha)
            new_data = ddpp.meanFieldExpansionSteadyState(order=order)
            test_data = data[key]
            assert approximately_equal(new_data, test_data) <= 1e-8

#generate_data()
#test_cache()
