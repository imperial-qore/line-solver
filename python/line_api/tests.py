import unittest
import linesolver as line
from . import *
from numpy import *

class LineSolverTests(unittest.TestCase):

    def setUp(self) :
        self.L = array([[10,5],[5,9]])
        self.N = array([100,100])
        self.Z = array([91,92])
    def test_pfqn_ca_1(self):
        [G, lG] = pfqn_ca(self.L,self.N,self.Z)
        self.assertEqual(lG, 540.9527921805867, 'Failed unit test')  # add assertion here

    def test_pfqn_aql(self):
        [XN,QN,UN,RN,AN,numIters] = pfqn_aql(self.L,self.N,self.Z)
        self.assertEqual(XN[0], 0.06159942822865749, 'Failed unit test')  # add assertion here

if __name__ == '__main__':
    unittest.main()

#%%
