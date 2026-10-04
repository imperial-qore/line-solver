"""
Basic tests to see if the derivatives are computed correctely
"""
import pickle
import sys
sys.path.append('..')
sys.path.append('.')
import src.rmf_tool as rmf
import numpy as np

import os
PWD=os.getcwd()
if PWD[-5:] == 'tests':
    CACHE_DIR = 'output_tests'
else:
    CACHE_DIR = 'tests/output_tests'


def sir_model():
    """
    this returns a density dependent population process of an SIR model
    """
    ddpp = rmf.DDPP()
    ddpp.add_transition([-1,1,0],lambda x:x[0]+2*x[0]*x[1])
    ddpp.add_transition([0,-1,+1],lambda x:x[1])
    ddpp.add_transition([1,0,-1],lambda x:3*x[2]**3)
    ddpp.set_initial_state([.3,.2,.5]) # We first need to define an initial stater
    return ddpp

def function(model, x0):
    """
    datas to be tested (from the model)
    """
    values = [model.defineDriftDerivativeQ(evaluate_at=x0),
              model.defineDriftSecondDerivativeQderivativesR(evaluate_at=x0)]
    return [array for mytuple in values for array in mytuple ]

def absolute_difference(new_data, old_data):
    """
    Takes two sequences of data and return the sum of the absolute differences between all
    """
    diff = 0
    assert len(new_data) == len(old_data)
    for new, old in zip(new_data, old_data):
        diff += np.sum(np.abs(new-old))
    return diff


def generate_data():
    """
    Generate a pickle file with the current version of the tool.
    """
    model = sir_model()
    data = dict([])
    for i in range(10):
        x0 = np.random.rand(3)
        x0 = x0/sum(x0)
        data[tuple(x0)] = function(model, x0)
        print(x0)
    with open('{}/drift_derivatives.pickle'.format(CACHE_DIR), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

def test_drift_derivatives():
    """
    Test if derivatives are correct
    """
    model = sir_model()
    with open('{}/drift_derivatives.pickle'.format(CACHE_DIR), 'rb') as f:
        data = pickle.load(f)
        for x0 in data:
            print(x0, 'OK')
            new_data = function(model, np.array(x0))
            test_data = data[x0]
            assert absolute_difference(new_data, test_data) <= 1e-8

#generate_data()
#test_drift_derivatives()
