import unittest
from line_solver import *
from numpy import *

class PFQNAPITests(unittest.TestCase):

    def setUp(self) :
        self.L = array([[10,5],[5,9]])
        self.N = array([100,100])
        self.Z = array([91,92])

    def test_pfqn_ca(self):
        [G, lG] = pfqn_ca(self.L,self.N,self.Z)
        self.assertEqual(lG, 549.1584415966641, 'Failed unit test')  # add assertion here

    def test_pfqn_bs(self):
        XN,QN,UN,RN,numIters = pfqn_bs(self.L,self.N,self.Z)
        self.assertEqual(XN[0], 0.06173674890340948, 'Failed unit test')  # add assertion here
        self.assertEqual(XN[1], 0.07587637803952278, 'Failed unit test')  # add assertion here
        self.assertEqual(QN[0][0], 72.416372121381727, 'Failed unit test')  # add assertion here
        self.assertEqual(QN[0][1], 44.606488689804834, 'Failed unit test')  # add assertion here
        self.assertEqual(QN[1][0], 21.965583728408010, 'Failed unit test')  # add assertion here
        self.assertEqual(QN[1][1], 48.412884530559083, 'Failed unit test')  # add assertion here
        self.assertEqual(UN[0][0], 0.6173674890340948, 'Failed unit test')  # add assertion here
        self.assertEqual(UN[0][1], 0.3793818901976139, 'Failed unit test')  # add assertion here
        self.assertEqual(UN[1][0], 0.3086837445170474, 'Failed unit test')  # add assertion here
        self.assertEqual(UN[1][1], 0.6828874023557050, 'Failed unit test')  # add assertion here
        self.assertEqual(RN[0][0], 1.1729864854834047E3, 'Failed unit test')  # add assertion here
        self.assertEqual(RN[0][1], 0.5878837372359818E3, 'Failed unit test')  # add assertion here
        self.assertEqual(RN[1][0], 0.35579430596150064E3, 'Failed unit test')  # add assertion here
        self.assertEqual(RN[1][1], 0.6380494928914713E3, 'Failed unit test')  # add assertion here

if __name__ == '__main__':
    unittest.main()

#%%
