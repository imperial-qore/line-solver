"""
Test that computes the refined mean field approximation for the two-choice model
(with order 1 and 2 and a few parameter)
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

def dChoiceModel(K, rho, d):
    ddpp = rmf.DDPP()

    # The vector 'e(i)' is a vector where the $i$th coordinate is equal to $1$ (the other being equal to $0$)
    def e(i):
        l = np.zeros(K)
        l[i] = 1
        return l

    # We then add the transitions : 
    for i in range(K):
        if i >= 1:
            ddpp.add_transition(e(i),eval('lambda x: {}*(x[{}]**{} - x[{}]**{} )'.format(rho, i-1, d, i, d)))
        if i < K-1:
            ddpp.add_transition(-e(i),eval('lambda x: (x[{}] - x[{}])'.format(i,i+1) ))
    ddpp.add_transition(e(0), lambda x : eval('{}*(1-x[0]**{})'.format(rho,d)))
    ddpp.add_transition(-e(K-1), lambda x : x[K-1])
    ddpp.set_initial_state(e(0))
    return ddpp

def generate_data():
    """
    Generate all data and store them in a pickle file
    (to be used one times when the test is initialized)
    """
    data = dict([])
    for rho in [0.6, 0.7, 0.8, 0.9]:
        for d in [2, 3]:
            for K in [5, 9, 15, 20]:
                for order in ([1, 2] if K <= 5 else [1]):
                    ddpp = dChoiceModel(K, rho, d)
                    data[(K, rho, d, order)] = ddpp.meanFieldExpansionSteadyState(order=order)
    with open('{}/d_choice.pickle'.format(CACHE_DIR), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

def test_two_choice():
    """
    Compare the new data with previously computed data.
    """
    with open('{}/d_choice.pickle'.format(CACHE_DIR), 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
        for key in data:
            (K,rho,d,order) = key
            print(key)
            ddpp = dChoiceModel(K, rho, d)
            new_data = ddpp.meanFieldExpansionSteadyState(order=order)
            test_data = data[key]
            assert approximately_equal(new_data, test_data) <= 1e-8

#generate_data()
#test_two_choice()
