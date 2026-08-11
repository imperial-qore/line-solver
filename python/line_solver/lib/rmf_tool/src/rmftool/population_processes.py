"""
This file provides the implementation of a class DDPP and a way to compute
the mean field and refined mean field approximation as defined in
https://hal.inria.fr/hal-01622054/document
and https://hal.inria.fr/hal-01891632/document

DDPP stands for Density Dependent Population Process
"""

import numpy as np
import random as rnd
import scipy.integrate as integrate
import sympy as sym
import scipy.linalg
import numpy.linalg
import matplotlib.pyplot as plt
import time as ti

import copy

from numpy import tensordot as tsdot
from sympy.utilities.lambdify import lambdify
from sympy import derive_by_array
import os
import sys

PATH = os.path.abspath(__file__)
DIR_PATH = os.path.dirname(PATH)
sys.path.append(DIR_PATH)

from refinedRefined_fixedPoint import computePiV, computePiVA
from refinedRefined_transientRegime import drift_r_vector, drift_rr_vector


### Additional imports for HetPP 
from het_refinedRefined_transientRegime import drift_r_vector
import scipy
import scipy.integrate as integrate
from scipy.interpolate import interp1d

# To compute the transient trajectories and fixed points



def eval_array_at(F, subs):
    """ evaluate an array of 'symbolic' expression using the dictionary 'subs'.

    Returns : an Array of float
    """
    myF = F.flatten()
    myF = np.array([float(f.evalf(subs=subs)) for f in myF])
    return myF.reshape(F.shape)


class RmfError(Exception):
    """Basic error class for this module
    """


class DimensionError(RmfError):
    pass


class NotImplemented(RmfError):
    pass


class NegativeRate(RmfError):
    pass


class InitialConditionNotDefined(RmfError):
    pass


class DDPP():
    """
    DDPP serves to define and study density depend population processes
    """

    def __init__(self, use_autograd=False):
        """
        Initialize an empty density dependent population process (without any transitions)

        Args:
            use_autograd (bool, optional): When True, the jacobian and hessian are computed
                by using autograd and not symbolic differentiation. This allows for more general 
                rate functions but is slower. Default=False
            param1 (str): Description of `param1`.
            param2 (:obj:`int`, optional): Description of `param2`. Multiple
                lines are supported.
        """
        self._list_of_transitions = []
        self._list_of_rate_functions = []
        self._x0 = None
        self._model_dimension = None
        self._use_autograd = use_autograd

    def add_transition(self, l, f):
        r"""
        Add a new transition of the form (\ell,\beta_\ell), where \beta_\ell(x)
        is the rate of the transition

        Args:
            l: A vector of changes (of size n)
            f: A function from R^n \to R, where f(x) is the rate at which the transition occurs

        Returns:
            True if successful, False otherwise

        Raises:
            DimensionError if l is not of the right size.
        """
        if self._model_dimension is not None and self._model_dimension != len(l):
            raise DimensionError
        self._model_dimension = len(l)
        self._list_of_transitions.append(np.array(l))
        self._list_of_rate_functions.append(f)

    def set_initial_state(self, x0):
        """
        Defines the initial state

        Args:
            x0 must be an array of the correct dimension
        """
        if self._model_dimension and self._model_dimension != len(x0):
            raise DimensionError
        self._model_dimension = len(x0)
        self._x0 = x0

    def simulate(self, N, time):
        """
        Simulates an realization of the stochastic process with N objects

        Returns:
            (T,X), where :
            - T is a 1-dimensional numpy array, where T[i] is the time of the i-th time step.
            - x is a 2-dimensional numpy array, where x[i,j] is the j-th coordinate of the system at time T[i]
        """
        if N == 'inf':
            return self.ode(time)
        if self._x0 is None:
            raise InitialConditionNotDefined
        nb_trans = len(self._list_of_transitions)
        t = 0

        # if fix!=-1:     seed(fix)
        x = np.array(self._x0)
        T = [0]
        X = [x]
        while t < time:
            L_rates = [self._list_of_rate_functions[i]
                       (x) for i in range(nb_trans)]
            if any(rate < -1e-14 for rate in L_rates):
                raise NegativeRate
            S = sum(L_rates)
            if S <= 1e-14:
                print('System stalled (total rate = 0)')
                t = time
            else:
                a = rnd.random()*S
                l = 0
                while a > L_rates[l]:
                    a -= L_rates[l]
                    l += 1

                x = x+(1./N)*self._list_of_transitions[l]

                t += rnd.expovariate(N*S)

            T.append(t)
            X.append(x)

        X = np.array(X)
        return(T, X)

    def ode(self, time, number_of_steps=1000):
        """Simulates the ODE (mean-field approximation)
        """
        if self._x0 is None:
            raise InitialConditionNotDefined

        def drift(x):
            return (sum([self._list_of_transitions[i]*self._list_of_rate_functions[i](x)
                         for i in range(len(self._list_of_transitions))], 0))

        T = np.linspace(0, time, number_of_steps)
        X = integrate.odeint(lambda x, t: drift(x), self._x0, T)

        return T, X

    def fixed_point(self):
        """
        Computes the fixed of the ODE (if this ODE has a fixed point starting from x0)
        """
        if self._x0 is None:
            print(
                'No initial condition given. We assume that the initial condition is "x0=[1,0,...]"')
            self._x0 = np.zeros(self._model_dimension)
            self._x0[0] = 1
        return self.ode(time=10000)[1][-1, :]

    def doTransitionsConserveSum(self):
        """
        This function tests if the transitions conserve the sum of the coordinates.

        Returns : True or False.
        """
        for l in self._list_of_transitions:
            if sum(l) != 0:
                return False
        return True

    def dimensionReduction(self, A):
        """
        """
        n = len(A)
        M = np.array([l for l in self._list_of_transitions])
        rank_of_transisions = np.linalg.matrix_rank(M)
        # If rank_of_transisions < n, this means that the stochastic process
        # evolves on a linear subspace of R^n.
        eigenvaluesOfJacobian = scipy.linalg.eig(A, left=False, right=False)
        rank_of_jacobian = np.linalg.matrix_rank(A)
        if sum(numpy.real(eigenvaluesOfJacobian) < 1e-8) < rank_of_transisions:
            # This means that there are less than "rank_of_transisions"
            # eigenvalues with <0 real part
            print("The Jacobian seems to be not Hurwitz")
        if rank_of_jacobian == n:
            return(np.eye(n), np.eye(n), n)
        U, _, _ = scipy.linalg.svd(A)
        return U.transpose(), U, rank_of_jacobian


    def test_for_linear_dependencies(self):
        """This function tests if there the transition conserve some subset of the space.

        Return (dim, subs) where
        * dim is the rank of the matrix of transitions
        * subs corresponds to the null space of the transitions"""

        M = np.array([l for l in self._list_of_transitions])
        u, d, v = np.linalg.svd(M, full_matrices=False)
        dim = sum(abs(d) > 1e-8)
        if dim < self._model_dimension:
            v = v[dim:]
            p, l, u = scipy.linalg.lu(v)
            subs = u
            for i in range(subs.shape[0]):
                j = np.argmax(abs(subs[i, :]) > 1e-5)
                subs[i, :] = subs[i, :] / subs[i, j]
        else:
            subs = np.array([])
        return(dim, subs)

    def _numericalJacobian(self, drift, fixedPoint, epsilon=1e-6):
        r""" Computes the Jacobian of the drift using a finite difference method.

        Returns : dF[i][j] = \partial f_i / \partial x_j evaluates in fixedPoint
        """
        dim = len(fixedPoint)

        def e(i):
            res = np.zeros(dim)
            res[i] = 1
            return res
        A = [(drift(fixedPoint+epsilon*e(j))-drift(fixedPoint -
                                                   epsilon*e(j)))/(2*epsilon) for j in range(dim)]
        return np.array(A).transpose()

    def _numericalHessian(self, drift, fixedPoint, epsilon=1e-4):
        r""" Computes the Jacobian of the drift using a finite difference method : 

        Returns : ddF[i][j][k] = \partial^2 f_i / (\partial x_j\partial x_k) evaluates in fixedPoint
        """
        dim = len(fixedPoint)

        def e(i):
            res = np.zeros(dim)
            res[i] = 1
            return(res)

        def ee(i, j):
            res = np.zeros(dim)
            res[i] = 1
            res[j] += 1
            return(res)
        ffixedPoint = drift(fixedPoint)
        ddB = [[drift(fixedPoint+epsilon*ee(i, j))
                for j in range(i+1)] for i in range(dim)]
        dB = [drift(fixedPoint+epsilon*e(i)) for i in range(dim)]
        B = [[(ddB[max(i, j)][min(i, j)] - dB[i] - dB[j] + ffixedPoint)/epsilon**2
              for i in range(dim)] for j in range(dim)]
        B = [[[B[j][k][l] for j in range(dim)]
              for k in range(dim)] for l in range(dim)]
        return(np.array(B))

    def defineDrift(self, evaluate_at=None):
        """Return a (lambdified) function F(.) that is the drift of the system

        Arg : evaluate_at=None : if not None, return the value of the drift at point "evaluate_at".
        """
        n = len(self._list_of_transitions[0])
        number_of_transitions = len(self._list_of_transitions)
        x = [sym.symbols('x[{}]'.format(i)) for i in range(n)]
        f = np.zeros(n)
        for l in range(number_of_transitions):
            f = f + self._list_of_transitions[l] * \
                self._list_of_rate_functions[l](x)
        F = sym.lambdify([x], [f[i] for i in range(n)])
        def computeF(x): return(np.array(F(x)))
        return computeF

    def defineDriftDerivativeQ_autograd(self, evaluate_at=None):
        """
        version of defineDriftDerivativeQ that uses 'JAX (successor of autograd)'
        """
        import jax.numpy as jnp
        import jax
        from jax.config import config
        config.update("jax_enable_x64", True) # we need to use double-float precision
        # (this could be used with single-precision but it is less accurate)

        n = len(self._list_of_transitions[0])
        number_of_transitions = len(self._list_of_transitions)
        def drift(x):
            f = jnp.zeros(n)
            for l in range(number_of_transitions):
                f = f + self._list_of_transitions[l] * \
                    self._list_of_rate_functions[l](x)
            return f
        jacobian = jax.jacobian(drift)
        hessian = jax.hessian(drift)
        def Q(x):
            myQ = np.zeros((n**2))
            for l in range(number_of_transitions):
                myQ += np.kron(self._list_of_transitions[l], self._list_of_transitions[l]
                               )*self._list_of_rate_functions[l](x)
            return myQ.reshape((n,n))
        if evaluate_at is not None:
            x = jnp.array(evaluate_at)
            return jacobian(x), hessian(x), Q(x)
        else:
            return jax.jit(jacobian), jax.jit(hessian), jax.jit(Q)


    def defineDriftDerivativeQ(self, evaluate_at=None):
        """Return three (lambdified) functions Fp(.), Fpp(.) and Q(.) that are the first two derivatives of the drift and the matrix Q. 

        Arg : evaluate_at=None : if not None, return the value of the functions evaluated at point "evaluate_at".
        """
        if self._use_autograd:
            return self.defineDriftDerivativeQ_autograd(evaluate_at=evaluate_at)
        
        n = len(self._list_of_transitions[0])
        number_of_transitions = len(self._list_of_transitions)
        x = [sym.symbols('x[{}]'.format(i)) for i in range(n)]
        f = np.zeros(n)
        for l in range(number_of_transitions):
            f = f + self._list_of_transitions[l] * \
                self._list_of_rate_functions[l](x)
        dF = np.array([[sym.diff(f[i], x[j]) for j in range(n)]
                       for i in range(n)])
        ddF = np.array([[[sym.diff(dF[i, j], x[k]) for j in range(n)]
                         for k in range(n)] for i in range(n)])

        if evaluate_at is not None:
            subs_dictionary = {x[i]: evaluate_at[i] for i in range(n)}
            Fp = eval_array_at(dF, subs_dictionary)
            Fpp = eval_array_at(ddF, subs_dictionary)
            Q = np.zeros((n**2))
            for l in range(number_of_transitions):
                Q = Q + \
                    np.kron(self._list_of_transitions[l], self._list_of_transitions[l]
                            )*self._list_of_rate_functions[l](evaluate_at)
            return(Fp, Fpp, Q.reshape((n, n)))

        q = np.zeros((n**2))
        for l in range(number_of_transitions):
            q = q + \
                np.kron(
                    self._list_of_transitions[l], self._list_of_transitions[l])*self._list_of_rate_functions[l](x)

        Fp = sym.lambdify([x], list(dF.reshape(n**2)))
        Fpp = sym.lambdify([x], list(ddF.reshape((n**3))))
        Q = sym.lambdify([x], [q[i] for i in range(n*n)])

        def computeFp(x): return(np.array(Fp(x)).reshape((n, n)))
        def computeFpp(x): return(np.array(Fpp(x)).reshape((n, n, n)))
        def computeQ(x): return(np.array(Q(x)).reshape((n, n)))
        return(computeFp, computeFpp, computeQ)

    def defineDriftSecondDerivativeQderivativesR(self, evaluate_at=None):
        """Return three (lambdified) functions Fppp(.), Fpppp(.) and Q(.) that are the 3rd and 4th derivatives of the drift, the first two derivatives of Q(.) and the matrix R(.)

        Arg : evaluate_at=None : if not None, return the value of the functions evaluated at point "evaluate_at".
        """
        n = len(self._list_of_transitions[0])
        number_of_transitions = len(self._list_of_transitions)
        x = [sym.symbols('x[{}]'.format(i)) for i in range(n)]
        f = np.zeros(n)
        for l in range(number_of_transitions):
            f = f + self._list_of_transitions[l] * \
                self._list_of_rate_functions[l](x)
        dddF = np.array([[sym.diff(f[i], x[j], x[k], x[l]) for j in range(n) for k in range(n) for l in range(n)]
                         for i in range(n)]).reshape((n**4))
        dddF_nonReshaped = dddF.reshape((n, n, n, n))
        ddddF = np.array([[sym.diff(dddF_nonReshaped[i, j, k, l], x[m]) for j in range(n) for k in range(n) for l in range(n) for m in range(n)]
                          for i in range(n)]).reshape((n**5))
        r = np.zeros((n**3))
        for l in range(number_of_transitions):
            r = r + np.kron(self._list_of_transitions[l], np.kron(
                self._list_of_transitions[l], self._list_of_transitions[l]))*self._list_of_rate_functions[l](x)
        q = np.zeros((n**2))
        for l in range(number_of_transitions):
            q = q + \
                np.kron(
                    self._list_of_transitions[l], self._list_of_transitions[l])*self._list_of_rate_functions[l](x)
        dQ = np.array([[[sym.diff(q[i], x[k]) for k in range(n)]
                        for i in range(n**2)]]).reshape(n**3)
        ddQ = np.array([[[sym.diff(q[i], x[k], x[l]) for k in range(n)
                          for l in range(n)] for i in range(n**2)]]).reshape(n**4)

        if evaluate_at is not None:
            subs_dictionary = {x[i]: evaluate_at[i] for i in range(n)}
            Fppp = eval_array_at(dddF.reshape((n, n, n, n)), subs_dictionary)
            Fpppp = eval_array_at(ddddF.reshape(
                (n, n, n, n, n)), subs_dictionary)
            Qp = eval_array_at(dQ.reshape((n, n, n)), subs_dictionary)
            Qpp = eval_array_at(ddQ.reshape((n, n, n, n)), subs_dictionary)
            R = eval_array_at(r.reshape((n, n, n)), subs_dictionary)
            return(Fppp, Fpppp, Qp, Qpp, R)

        Fppp = sym.lambdify([x], [dddF[i] for i in range(n**4)])
        Fpppp = sym.lambdify([x], [ddddF[i] for i in range(n**5)])
        Qp = sym.lambdify([x], [dQ[i] for i in range(n**3)])
        Qpp = sym.lambdify([x], [ddQ[i] for i in range(n**4)])
        R = sym.lambdify([x], [r[i] for i in range(n**3)])
        def computeFppp(x): return(np.array(Fppp(x)).reshape((n, n, n, n)))
        def computeFpppp(x): return(
            np.array(Fpppp(x)).reshape((n, n, n, n, n)))

        def computeQp(x): return(np.array(Qp(x)).reshape((n, n, n)))
        def computeQpp(x): return(np.array(Qpp(x)).reshape((n, n, n, n)))
        def computeR(x): return(np.array(R(x)).reshape((n, n, n)))
        return(computeFppp, computeFpppp, computeQp, computeQpp, computeR)

    def meanFieldExpansionTransient(self, order=1, time=10):
        """ Computes the transient values of the mean field approximation or its O(1/N^{order})-expansions

        Args:
           - order : can be 0 (mean field approx.), 1 (O(1/N)-expansion) or 2 (O(1/N^2)-expansion)

        Returns : (T,XVW) or (T,XVWABCD), where T is a time interval and XVW is a (2d+d^2)*number_of_steps matrix (or XVWABCD is a (3n+2n^2+n^3+n^4) x number_of_steps matrix), where : 
        * XVW[0:n,:]                 is the solution of the ODE (= mean field approximation)
        * XVW[n:2*n,:]               is V(t) (= 1st order correction)
        * XVW[2*n:2*n+n**2,:]        is W(t)
        * XVWABCD[2*n+n**2,3*n+n**2] is A(t) (= the 2nd order correction)
        """
        n = len(self._list_of_transitions[0])
        t_start = ti.time()

        # We first defines the function that will be used to compute the drift (using symbolic computation)
        computeF = self.defineDrift()
        if (order >= 1):  # We need 2 derivatives and Q to get the O(1/N)-term
            computeFp, computeFpp, computeQ = self.defineDriftDerivativeQ()
        if (order >= 2):  # We need the next 2 derivatives of F and Q + the tensor R
            computeFppp, computeFpppp, computeQp, computeQpp, computeR = self.defineDriftSecondDerivativeQderivativesR()
        print('time to compute drift=', ti.time()-t_start)

        if order == 0:
            X_0 = self._x0
            Tmax = time
            T = np.linspace(0, Tmax, 1000)
            numericalInteg = integrate.solve_ivp(lambda t, x: computeF(x), [
                                                 0, Tmax], X_0, t_eval=T, rtol=1e-6)
            return(numericalInteg.t, numericalInteg.y.transpose())
        if order == 1:
            XVW_0 = np.zeros(2*n+n**2)
            XVW_0[0:n] = self._x0

            Tmax = time
            T = np.linspace(0, Tmax, 1000)

            numericalInteg = integrate.solve_ivp(lambda t, x:
                                                 drift_r_vector(
                                                     x, n, computeF, computeFp, computeFpp, computeQ),
                                                 [0, Tmax], XVW_0, t_eval=T, rtol=1e-6)
            XVW = numericalInteg.y.transpose()
            X = XVW[:, 0:n]
            V = XVW[:, n:2*n]
            return(numericalInteg.t, X, V, XVW)
        elif order == 2:
            XVWABCD_0 = np.zeros(3*n+2*n**2+n**3+n**4)
            XVWABCD_0[0:n] = self._x0

            Tmax = time
            T = np.linspace(0, Tmax, 1000)

            numericalInteg = integrate.solve_ivp(lambda t, x:
                                                 drift_rr_vector(x, n, computeF, computeFp, computeFpp, computeQ,
                                                                 computeFppp, computeFpppp, computeQp, computeQpp, computeR),
                                                 [0, Tmax], XVWABCD_0, t_eval=T, rtol=1e-6)
            XVWABCD = numericalInteg.y.transpose()
            X = XVWABCD[:, 0:n]
            V = XVWABCD[:, n:2*n]
            A = XVWABCD[:, 2*n+n**2:3*n+n**2]
            return(numericalInteg.t, X, V, A, XVWABCD)
        else:
            print(
                "order must be 0 (mean field), 1 (refined of order O(1/N)) or 2 (refined order 1/N^2)")

    def reduceDimensionFpFppQ(self, Fp, Fpp, Q):
        P, Pinv, rank = self.dimensionReduction(Fp)
        Fp = (P@Fp@Pinv)[0:rank, 0:rank]
        Fpp = tsdot(tsdot(tsdot(P, Fpp, 1), Pinv, 1), Pinv,
                    axes=[[1], [0]])[0:rank, 0:rank, 0:rank]
        Q = (P@Q@P.transpose())[0:rank, 0:rank]
        return(Fp, Fpp, Q, P, Pinv, rank)

    def expandDimensionVW(self, V, W, Pinv):
        rank = len(V)
        return(Pinv[:, 0:rank]@ V, Pinv[:, 0:rank]@W@Pinv.transpose()[0:rank, :])

    def reduceDimensionFppFppppQpQppR(self, Fppp, Fpppp, Qp, Qpp, R, P, Pinv, rank):
        FpppOLD = Fppp
        FppppOLD = Fpppp
        QpOLD = Qp
        Fppp = tsdot(tsdot(tsdot(tsdot(P, Fppp, 1), Pinv, 1), Pinv, axes=[
                     [1], [0]]), Pinv, axes=[[2], [0]])[0:rank, 0:rank, 0:rank, 0:rank]
        Fpppp = tsdot(tsdot(tsdot(tsdot(tsdot(P, Fpppp, 1), Pinv, 1), Pinv, axes=[[1], [0]]), Pinv, axes=[
                      [2], [0]]), Pinv, axes=[[3], [0]])[0:rank, 0:rank, 0:rank, 0:rank, 0:rank]
        Qp_new = np.transpose(tsdot(tsdot(P, np.transpose(
            Qp, axes=[0, 2, 1]), 1), P.transpose(), 1), axes=[0, 2, 1])
        Qp_final = tsdot(Qp_new, Pinv, 1)[0:rank, 0:rank, 0:rank]
        Qpp_new = np.transpose(tsdot(tsdot(P, np.transpose(
            Qpp, axes=[0, 3, 2, 1]), 1), P.transpose(), 1), axes=[0, 3, 2, 1])
        Qpp = tsdot(tsdot(Qpp_new, Pinv, axes=[[2], [0]]), Pinv, axes=[
                    [3], [0]])[0:rank, 0:rank, 0:rank, 0:rank]
        R = tsdot(tsdot(tsdot(P, R, 1), P.transpose(), axes=[
                  [1], [0]]), P.transpose(), 1)[0:rank, 0:rank, 0:rank]
        return(Fppp, Fpppp, Qp_final, Qpp, R)

    def expandDimensionABCD(self, A, B, C, D, Pinv):
        rank = len(A)
        Pinv = Pinv[:, 0:rank]
        return(Pinv@A, Pinv@B@Pinv.transpose(),
               tsdot(tsdot(tsdot(C, Pinv, axes=[[2], [1]]), Pinv, axes=[
                     [1], [1]]), Pinv, axes=[[0], [1]]),
               tsdot(tsdot(tsdot(tsdot(D, Pinv, axes=[[3], [1]]), Pinv, axes=[[2], [1]]), Pinv, axes=[[1], [1]]), Pinv, axes=[[0], [1]]))

    def meanFieldExpansionSteadyState(self, order=1):
        """This code computes the O(1/N) and O(1/N^2) expansion of the mean field approximaiton
        (the term "V" is the "V" of Theorem~1 of https://hal.inria.fr/hal-01622054/document.

        Note : Probably less robust and slower that theoretical_V

        """
        pi = self.fixed_point()
        if order == 0:
            return pi
        if (order >= 1):  # We need 2 derivatives and Q to get the O(1/N)-term
            Fp, Fpp, Q = self.defineDriftDerivativeQ(evaluate_at=pi)
            Fp, Fpp, Q, P, Pinv, rank = self.reduceDimensionFpFppQ(Fp, Fpp, Q)
            if order == 1:
                pi, V, (V, W) = computePiV(pi, Fp, Fpp, Q)
                V, W = self.expandDimensionVW(V, W, Pinv)
                return(pi, V, (V, W))
        if (order >= 2):  # We need the next 2 derivatives of F and Q + the tensor R
            if len(self._x0) >= 11:
                print("*Warning* The computation time grows quickly with the number of dimensions",
                      "(probably around {:1.0f}sec. for this model)".format(0.0001*len(self._x0)**5))
            Fppp, Fpppp, Qp, Qpp, R = self.defineDriftSecondDerivativeQderivativesR(
                evaluate_at=pi)
            Fppp, Fpppp, Qp, Qpp, R = self.reduceDimensionFppFppppQpQppR(
                Fppp, Fpppp, Qp, Qpp, R, P, Pinv, rank)
            pi, V, A, (V, W, A, B, C, D) = computePiVA(
                pi, Fp, Fpp, Fppp, Fpppp, Q, Qp, Qpp, R)
            V, W = self.expandDimensionVW(V, W, Pinv)
            A, B, C, D = self.expandDimensionABCD(A, B, C, D, Pinv)
            if order == 2:
                return(pi, V, A, (V, W, A, B, C, D))
            else:
                print('order should be 0, 1 or 2')

    def theoretical_V(self, symbolic_differentiation=True):
        """This code computes the constant "V" of Theorem~1 of https://hal.inria.fr/hal-01622054/document

        Note : for now this function does not support rates that depend on N (i.e. C=0)

        This function is deprecated (meanFieldExpansionSteadyState(order=1) should be called instead)

        Args:
            symbolic_differentiation (bool, default=True) : when True, the derivative are computed using Sympy. When false, they are computed using a finite difference method.
        """
        n = self._model_dimension
        number_transitions = len(self._list_of_transitions)
        fixedPoint = self.fixed_point()

        Var = np.array([sym.symbols('x_{}'.format(i)) for i in range(n)])
        f_x = np.zeros(n)
        for i in range(number_transitions):
            f_x = f_x + \
                self._list_of_transitions[i] * \
                self._list_of_rate_functions[i](Var)

        # The following code attempts to reduce the number of dimensions by testing
        # if the transitions conserve some linear combinations of the coordinates
        (dim, subs) = self.test_for_linear_dependencies()
        variables = [i for i in range(n)]
        for sub in range(subs.shape[0]):
            j = np.argmax(abs(subs[sub, :]) > 1e-5)
            variables.remove(j)
            for i in range(n):
                f_x[i] = f_x[i].subs(Var[j], sum([self._x0[k]*subs[sub, k] for k in range(j, n)])
                                     - sum(np.array([Var[k]*subs[sub, k] for k in range(j+1, n)])))
        Var = [Var[i] for i in variables]

        if (symbolic_differentiation):
            dF = [[sym.diff(f_x[variables[i]], Var[j])
                   for j in range(dim)] for i in range(dim)]
            subs_dictionary = {Var[i]: fixedPoint[variables[i]]
                               for i in range(dim)}
            A = np.array([[float(dF[i][j].evalf(subs=subs_dictionary))
                           for j in range(dim)]
                          for i in range(dim)])
            B = np.array([[[float(sym.diff(dF[j][k], Var[l]).evalf(subs=subs_dictionary))
                            for l in range(dim)]
                           for k in range(dim)]
                          for j in range(dim)])
        else:
            drift = lambdify([Var], [f_x[variables[i]] for i in range(dim)])
            def drift_array(x): return np.array(drift(x))
            fixedPointProj = np.array(
                [fixedPoint[variables[i]] for i in range(dim)])
            A = self._numericalJacobian(drift_array, fixedPointProj)
            B = self._numericalHessian(drift_array, fixedPointProj)

        Q = np.zeros((dim, dim))
        for l in range(number_transitions):
            v = [self._list_of_transitions[l][variables[p]]
                 for p in range(dim)]
            Q += np.kron(v, v).reshape(dim, dim) * \
                self._list_of_rate_functions[l](fixedPoint)

        W = scipy.linalg.solve_continuous_lyapunov(A, Q)
        V = tsdot(scipy.linalg.inv(A), tsdot(B, W/2, 2), 1)

        # We now attemps to reconstruct the full C if the number of dimensions was reduced.
        if dim < n:
            newV = np.zeros(n)
            for i in range(dim):
                newV[variables[i]] = V[i]
            for sub in reversed(range(len(subs))):
                j = np.argmax(abs(subs[sub, :]) > 1e-5)
                newV[j] = -sum(np.array([newV[k]*subs[sub, k]
                                         for k in range(j+1, n)]))
            V = newV
        return(np.array(V))

    def _batch_meanConfidenceInterval(self, T, X):
        n = len(T)
        if len(T) < 300:  # We do not do batch if we have less than 200 samples
            return(np.sum(np.diff(T[int(n/2):n])*X[int(n/2):n-1]) / (T[n-1]-T[int(n/2)]), None)
        # Otherwise we split the data into 100 batchs of size
        Y = [np.sum(np.diff(T[int(2*i*n/300):int((2*i+1)*n/300)]) *
                    X[int(2*i*n/300):int((2*i+1)*n/300)-1] /
                    (T[int((2*i+1)*n/300)]-T[int(2*i*n/300)])) for i in range(50, 150)]
        return (np.mean(Y), 2*np.std(Y)/10)

    def steady_state_simulation(self, N, time=1000):
        """Generates a sample of E[X] in steady-state.

        The expectation is computed by performing one simulation from
        t=0 to t=time and averaging the values over time/2 .. time
        (the implicit assumption is that the system should be roughly
        in steady-state at time/2).

        """
        T, X = self.simulate(N, time)
        n = len(T)
        n2 = int(n/2)
        result = np.array([self._batch_meanConfidenceInterval(
            T, X[:, i]) for i in range(len(X[0, :]))])
        return(result[:, 0], result[:, 1])

    def compare_refinedMF(self, N, time=1000):
        """Compare E[X] with its mean-field and refined mean-field approximation

        Args : 
             N (int) : system's size
             time : Computes the expectation by

        Return : (Xm,Xrmf,Xs,Vs) where: 
             Xm is the fixed point of the mean-field approximation
             Xrmf = Xm + V/N
             Xs is an approximation E[X] (computed by simulation)
             Vs is an estimation of the variance
        """
        Xm = self.fixed_point()
        V = self.theoretical_V()
        Xrmf = Xm+V/N
        Xs, Vs = self.steady_state_simulation(N, time)
        return(Xm, Xrmf, Xs, Vs)

    def plot_ODE_vs_simulation(self, N, time=10):
        """Plot the ODE and the simulation on the same graph

        The system is simulated from t=0 to t=time (both for the
        stochastic and ODE version), starting from x0.

        Args: 
             N (int) : system's size
             time : the time until the simulation should be run

        """
        T, X = self.simulate(N, time=10)
        plt.plot(T, X)
        plt.gca().set_prop_cycle(None)
        T, X = self.ode(time=10)
        plt.plot(T, X, '--')
        plt.legend(['x_{}'.format(i) for i in range(self._model_dimension)])
        plt.xlabel('Time')
        plt.ylabel('x_{i}')
        plt.show()


class HomPP(DDPP):
    def __init__(self):
        super().__init__()

    def add_rate_tensors(self, unilateral_transition, pairwise_transitions):
        r"""
        Rate tensors must have the form:
        A[s, s_prime] specifying the unilateral transitions of objects from s to s_prime
        B[s, s1, s_prime, s1_prime] pairwise transitions of objects from  s to s_prime and s1 to s1_prime

        A has dimensions S x S
        B has dimensions S x S x S x S
        """

        self.unilateral_transition = unilateral_transition
        self.pairwise_transitions = pairwise_transitions

        trans_indices_unilateral = np.nonzero(self.unilateral_transition)
        nr_unilateral_trans = trans_indices_unilateral[0].shape[0]

        # add unilateral transitions
        for i in range(nr_unilateral_trans):
            s, s_prime = trans_indices_unilateral[0][i], trans_indices_unilateral[1][i]
            _transition_vector = - self.e(s) + self.e(s_prime)
            # _transition_rate = lambda x: self.unilateral_transition[s, s_prime] * x[s]
            _transition_rate = eval(
                'lambda x: {} * x[{}]'.format(self.unilateral_transition[s, s_prime], s)
            )
            self.add_transition(_transition_vector, _transition_rate)

        trans_indices_pairwise = np.nonzero(self.pairwise_transitions)
        nr_pair_trans = trans_indices_pairwise[0].shape[0]

        # add pairwise transitions
        for i in range(nr_pair_trans):
            s, s1 = trans_indices_pairwise[0][i], trans_indices_pairwise[1][i]
            s_prime, s1_prime = trans_indices_pairwise[2][i], trans_indices_pairwise[3][i]
            _transition_vector = - self.e(s) + self.e(s_prime) - self.e(s1) + self.e(s1_prime)
            _transition_rate = eval('lambda x: {} * x[{}] * x[{}]'.format(
                self.pairwise_transitions[s, s1, s_prime, s1_prime], s, s1)
            )
            # _transition_rate = lambda x: self.pairwise_transitions[s, s1, s_prime, s1_prime] * x[s] * x[s1]
            self.add_transition(_transition_vector, copy.deepcopy(_transition_rate))

    def e(self, s):
        _e = np.zeros(shape=self.unilateral_transition.shape[0])
        _e[s] = 1
        return _e


class HetPP(DDPP):
    def __init__(self):
        # initialize the heterogeneous model
        super().__init__()
        self._model_dimension_N = None
        self._model_dimension_S = None

        self.ddF = None

    def add_rate_tensors(self, A, B):
        r"""
        Rate tensors must have the form:
        A[k, s, s_prime] specifying the unilateral transition of object k from s to s_prime
        B[k, k1, s, s1, s_prime, s1_prime] pairwise transition of object k from s to s_prime and k1 from s1 to s1_prime

        A has dimensions N x S x S
        B has dimensions N x N x S x S x S x S
        """

        # set model dimensions
        if self._model_dimension_N is not None and self._model_dimension_N != A.shape[0]:
            if self._model_dimension_S is not None and self._model_dimension_S != B.shape[1]:
                raise DimensionError

        # infer model dimensions
        self._model_dimension_N = A.shape[0]
        self._model_dimension_S = A.shape[1]

        self._model_dimension = (self._model_dimension_N, self._model_dimension_S)

        # save rate tensors in model
        if np.sum(A < -1e-14) >= 1:
            raise NegativeRate("Rate tensor A has negative rates.")
        if np.sum(B < -1e-14) >= 1:
            raise NegativeRate("Rate tensor B has negative rates")
        self.A = A
        self.B = B

    def simulate(self, time, seed=None):
        """Simulates a sample trajectory of the system.
        Arguments:
            time - specifies up to which time the system should be simulated
            seed - sets the seed for the random generator; standard value - None
            
        Returns:
            T - jump times of the process
            X - values of the process
        """
        if self._x0 is None:
            raise InitialConditionNotDefined

        trans_indices_A = np.nonzero(self.A)
        nr_trans_A = trans_indices_A[0].shape[0]

        trans_indices_B = np.nonzero(self.B)
        nr_trans_B = trans_indices_B[0].shape[0]

        x = copy.copy(self._x0)
        t = 0
        T = [0]
        X = [copy.copy(x)]

        def transition_rates(_x):

            _trans_rates_A = np.zeros(shape=trans_indices_A[0].shape)

            for j in range(nr_trans_A):
                _k, _s, _s_prime = trans_indices_A[0][j], trans_indices_A[1][j], trans_indices_A[2][j]
                _trans_rates_A[j] = self.A[_k, _s, _s_prime] * _x[_k, _s]

            _trans_rates_B = np.zeros(shape=trans_indices_B[0].shape)

            for j in range(nr_trans_B):
                _k, _k1 = trans_indices_B[0][j], trans_indices_B[1][j]
                _s, _s1 = trans_indices_B[2][j], trans_indices_B[3][j]
                _s_prime, _s1_prime = trans_indices_B[4][j], trans_indices_B[5][j]
                _trans_rates_B[j] = self.B[_k, _k1, _s, _s1, _s_prime, _s1_prime] * _x[_k, _s] * _x[_k1, _s1]

            return _trans_rates_A, _trans_rates_B

        # set random generator to seed
        rnd.seed(seed)

        while t < time:
            # transition rates evaluated at x
            trans_rates_A, trans_rates_B = transition_rates(x)
            S = np.sum(trans_rates_A) + np.sum(trans_rates_B)

            if S <= 1e-12:
                print('System stalled (total rate = 0)')
                t = time
            else:
                a = rnd.random() * S
                # selected transition
                transition_selected = False
                for i in range(trans_rates_A.shape[0]):
                    if a > trans_rates_A[i]:
                        a -= trans_rates_A[i]
                    else:
                        transition_selected = True
                        k, s, s_prime = trans_indices_A[0][i], trans_indices_A[1][i], trans_indices_A[2][i]
                        # add unilateral transition to state vector
                        x += (- self.e(k, s) + self.e(k, s_prime))
                        break

                if transition_selected is False:
                    for i in range(trans_rates_B.shape[0]):
                        if a > trans_rates_B[i]:
                            a -= trans_rates_B[i]
                        else:
                            k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
                            s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
                            s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
                            # add unilateral transition to state vector
                            x += (- self.e(k, s) + self.e(k, s_prime) - self.e(k1, s1) + self.e(k1, s1_prime))
                            break

                t += rnd.expovariate(S)

            T.append(t)
            X.append(copy.copy(x))

        X = np.array(X)
        return T, X

    def sampleMean(self, time=10, samples=1000, time_steps=2000, seeds=None):
        """
        Calculates the sample mean for a specified number of sample trajectories of the system.
        The sample mean is derived from interpolated sample trajectories which are evaluated in
        uniformly spaced time steps up to a time t.

        Returns:
            interpolation_times - array of times at which the sample mean and squared error were evaluated
            sample_mean - array with sample mean values
        """
        # initialize the mean values (simulation data come as array nr_items*nr_caches)
        sample_mean = np.zeros(shape=(time_steps, self._model_dimension_N, self._model_dimension_S))

        # print("Calculating mean for " + str(nr_simulations) + " simulations.")
        interpolation_times = np.linspace(0, time, time_steps)

        print("Mean Calculation; Nr. Simulations done:")
        for i in range(samples):
            if ((i + 1) % 10) == 0:
                print(i + 1, ' ', end='')
            if seeds is None:
                _T, _X = self.simulate(time=time, seed=None)
            else:
                _T, _X = self.simulate(time=time, seed=seeds[i])
            # append new interpolation to list important that result from interpolations need to be transposed again
            # to obtain (nr_item, state) sized array
            interpolation = interp1d(_T, _X, axis=0)
            for k, time in enumerate(interpolation_times):
                sample_mean[k] += interpolation(time)
        sample_mean *= 1 / float(samples)

        return interpolation_times, sample_mean

    def sampleMeanVariance(self, time=10, samples=1000, time_steps=2000, seeds=None):
        """
        Calculates the sample mean and squared error for a specified number of sample trajectories of the system.
        Both, sample mean and squared error are derived from interpolated sample trajectories which are evaluated in
        uniformly spaced time steps up to a time t.

        Returns:
            interpolation_times - array of times at which the sample mean and squared error were evaluated
            sample_mean - array with sample mean values
            sample_squared_error - array with sample squared errors
        """

        # initialize the mean and squared error values (simulation data is returned as array nr_items*nr_caches)
        sample_mean = np.zeros(shape=(time_steps, self._model_dimension_N, self._model_dimension_S))
        sample_squared_error = np.zeros(shape=(time_steps, self._model_dimension_N, self._model_dimension_S))

        # specify times for interpolation
        interpolation_times = np.linspace(0, time, time_steps)

        print("Mean Calculation; Nr. Simulations done:")
        for i in range(samples):
            if ((i + 1) % 10) == 0:
                print(i + 1, ' ', end='')
            if seeds is None:
                _T, _X = self.simulate(time=time, seed=None)
            else:
                _T, _X = self.simulate(time=time, seed=seeds[i])
            # append new interpolation to list important that result from interpolations need to be transposed again
            # to obtain (nr_item, state) sized array
            interpolation = interp1d(_T, _X, axis=0)
            for k, time in enumerate(interpolation_times):
                sample_mean[k] += interpolation(time)
        sample_mean *= 1 / float(samples)

        print("\nVariance; Nr. Simulations done:")
        for i in range(samples):
            if ((i + 1) % 10) == 0:
                print(i + 1, ' ', end='')
            if seeds is None:
                _T, _X = self.simulate(time=time, seed=None)
            else:
                _T, _X = self.simulate(time=time, seed=seeds[i])
            # append new interpolation to list important that result from interpolations need to be transposed again
            # to obtain (nr_item, state) sized array
            interpolation = interp1d(_T, _X, axis=0)
            for k, time in enumerate(interpolation_times):
                sample_squared_error[k] += np.power((interpolation(time) - sample_mean[k]), 2)
        sample_squared_error *= 1 / float(samples - 1)

        return interpolation_times, sample_mean, sample_squared_error

    def defineDrift(self, evaluate_at=None):
        """
        Defines the drift vector given the transition tensors A and B.
        The drift is of the dimension N x S .

        Returns:
            drift - drift of the system at a given point ( numpy array of size N x S )
        """

        # raise error if evaluate_at is None
        if evaluate_at is None:
            raise NotImplemented("A lambdified version of the drift is currently only available for the homogeneous "
                                 "implementation.")

        # drift vector
        drift = np.zeros(shape=(self._model_dimension_N, self._model_dimension_S))

        x = evaluate_at

        # use non zero rate entries to filter actual transitions
        trans_indices_A = np.nonzero(self.A)
        nr_pair_trans_A = trans_indices_A[0].shape[0]

        for i in range(nr_pair_trans_A):
            k, s, s_prime = trans_indices_A[0][i], trans_indices_A[1][i], trans_indices_A[2][i]
            drift += self.A[k, s, s_prime] * x[k, s] * (- self.e(k, s) + self.e(k, s_prime))

        trans_indices_B = np.nonzero(self.B)
        nr_pair_trans_B = trans_indices_B[0].shape[0]

        for i in range(nr_pair_trans_B):
            k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
            s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
            s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
            drift += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s] * x[k1, s1] * \
                     (- self.e(k, s) + self.e(k, s_prime) - self.e(k1, s1) + self.e(k1, s1_prime))

        return drift

    def set_initial_state(self, x0):
        r"""
        Sets the initial state of the system to x0.
        """
        # make sure initial state dimension match with from tensor inferred dimensions
        if self._model_dimension_N is not None and self._model_dimension_S is not None:
            if self._model_dimension_N != x0.shape[0] or self._model_dimension_S != x0.shape[1]:
                raise DimensionError("The model dimensions do not match. Make sure the model is properly initialized.")
        if self._model_dimension_N is None and self._model_dimension_S is None:
            self._model_dimension_N, self._model_dimension_S = x0.shape[0], x0.shape[1]
        self._x0 = x0

    def defineDriftDerivativeQ_autograd(self, evaluate_at=None):
        """
        Autograd and Jax methods are not available for the heteroegeneous model.
        """
        raise NotImplementedError("The defineDriftDerivativeQ_autograd method is only available for "
                                  "the homogeneous model.")

    def defineDriftDerivativeQ(self, evaluate_at=None):
        """
        Defines the drift derivatives and Q matrix at a given point. Used to calculate the refinement term.
        """

        if evaluate_at is None:
            raise NotImplemented("A lambdified version of the drift is currently only "
                                 "available for the homogeneous implementations.")

        # redefine evaluation point to define the derivatives and Q more comprehensively
        x = evaluate_at

        # use non zero rate entries to filter actual transitions
        trans_indices_A = np.nonzero(self.A)
        nr_unit_trans_A = trans_indices_A[0].shape[0]

        trans_indices_B = np.nonzero(self.B)
        nr_pair_trans_B = trans_indices_B[0].shape[0]

        # initialize Q matrix with dimension (N x S x N x S)
        Q = np.zeros(shape=(self._model_dimension_N, self._model_dimension_S,
                            self._model_dimension_N, self._model_dimension_S))

        # add transitions to the Q tensor
        # unilateral transitions
        for i in range(nr_unit_trans_A):
            # transition indices
            k, s, s_prime = trans_indices_A[0][i], trans_indices_A[1][i], trans_indices_A[2][i]
            # transition tensor
            transition_step = - self.e(k, s) + self.e(k, s_prime)
            transition_tensor = np.tensordot(transition_step, transition_step, axes=0)
            # add weighted transition to Q
            Q += self.A[k, s, s_prime] * x[k, s] * transition_tensor

        # pairwise transitions
        for i in range(nr_pair_trans_B):
            # transition indices
            k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
            s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
            s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
            # transition tensor
            transition_step = (- self.e(k, s) + self.e(k, s_prime) - self.e(k1, s1) + self.e(k1, s1_prime))
            transition_tensor = np.tensordot(transition_step, transition_step, axes=0)
            # add weighted transition to Q
            Q += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s] * x[k1, s1] * transition_tensor

        # define first derivative of drift
        # Jakobian has the dimensions (N x S) x (N x S)
        dF = np.zeros(shape=(self._model_dimension_N, self._model_dimension_S,
                             self._model_dimension_N, self._model_dimension_S))
        # unilateral transitions
        for i in range(nr_unit_trans_A):
            # transition indices
            k, s, s_prime = trans_indices_A[0][i], trans_indices_A[1][i], trans_indices_A[2][i]
            dF[k, s, k, s] -= self.A[k, s, s_prime]
            dF[k, s_prime, k, s] += self.A[k, s, s_prime]

        # pairwise transitions
        for i in range(nr_pair_trans_B):
            k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
            s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
            s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
            # define derivatives
            dF[k, s, k, s] -= self.B[k, k1, s, s1, s_prime, s1_prime] * x[k1, s1]
            dF[k1, s1, k, s] -= self.B[k, k1, s, s1, s_prime, s1_prime] * x[k1, s1]
            dF[k, s, k1, s1] -= self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s]
            dF[k1, s1, k1, s1] -= self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s]

            dF[k, s_prime, k, s] += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k1, s1]
            dF[k1, s1_prime, k, s] += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k1, s1]
            dF[k, s_prime, k1, s1] += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s]
            dF[k1, s1_prime, k1, s1] += self.B[k, k1, s, s1, s_prime, s1_prime] * x[k, s]

        # define first derivative of drift
        # Hessian has the dimensions (N x S) x (N x S) x (N x S)

        if self.ddF is None:

            ddF = np.zeros(shape=(self._model_dimension_N, self._model_dimension_S,
                                  self._model_dimension_N, self._model_dimension_S,
                                  self._model_dimension_N, self._model_dimension_S))
            # pairwise transitions ( no unilateral transitions)
            for i in range(nr_pair_trans_B):
                k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
                s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
                s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
                # define derivatives
                # drift / dx[k,s] dx[k1,s1]
                ddF[k, s, k, s, k1, s1] -= self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k, s_prime, k, s, k1, s1] += self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k1, s1, k, s, k1, s1] -= self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k1, s1_prime, k, s, k1, s1] += self.B[k, k1, s, s1, s_prime, s1_prime]
                # drift / dx[k1,s1] dx[k,s] - permutation of derivative direction
                ddF[k, s, k1, s1, k, s] -= self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k, s_prime, k1, s1, k, s] += self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k1, s1, k1, s1, k, s] -= self.B[k, k1, s, s1, s_prime, s1_prime]
                ddF[k1, s1_prime, k1, s1, k, s] += self.B[k, k1, s, s1, s_prime, s1_prime]

            self.ddF = ddF
        else:
            ddF = self.ddF

        return dF, ddF, Q

    def meanFieldExpansionTransient(self, order=1, time=10):
        """
        Computes the transient values of the mean field approximation or its O(1/N^{order})-expansions

        Args:
           - order : can be 0 (mean field approx.) or 1 (O(1/N)-expansion);


        Returns :
        if order = 0 
            T, X  - where T is a time interval, X is the solution to mean field approximation
        if order = 1
           T, X, V, XVW - where T is a time interval, X is the solution to the mean field approximation,
            V is the solution to the refinement term and XVW is a (2(N*S)+(N*S)^2)*number_of_steps matrix  where :
                * XVW[0:n,:]                 is the solution of the ODE (= mean field approximation, equal to X)
                * XVW[n:2*n,:]               is V(t) (= 1st order correction, equal to V)
                * XVW[2*n:2*n+n**2,:]        is W(t)
        """

        if order >= 2:
            raise NotImplementedError("Second order methods can not be implemented for polynomial drift of order two.")

        # set model dimensions
        N = self._model_dimension_N
        S = self._model_dimension_S

        if order == 0:
            # return mean field model
            T, X = self.ode(time=time, number_of_steps=1000)
            return T, X

        if order == 1:
            # allocate N*S values for the mean field X, N*S values for the expansion term V
            # and (N*S)**2 values for the expansion term W
            XVW_0 = np.zeros(2 * N * S + (N * S) ** 2)
            XVW_0[0:N * S] = self._x0.flatten()

            # max time of the integration
            t_max = time
            # specify time steps
            T = np.linspace(0, t_max, 1000)

            # calculating the solution to the ivp
            numerical_integration = integrate.solve_ivp(lambda t, x:
                                                        drift_r_vector(
                                                            x, N, S, self.defineDrift, self.defineDriftDerivativeQ
                                                        ),
                                                        [0, t_max], XVW_0, t_eval=T, rtol=1e-8)

            # separating and reshaping X and V values
            XVW = numerical_integration.y.transpose()
            T = numerical_integration.t

            X = XVW[:, 0:(N * S)]
            # reshape X values into (N x S) format
            X = np.array([X[i, :].reshape((N, S)) for i in range(X.shape[0])])

            V = XVW[:, (N * S):2 * (N * S)]
            # reshape V values into (N x S) format
            V = np.array([V[i, :].reshape((N, S)) for i in range(V.shape[0])])

            return T, X, V, XVW

    def meanFieldExpansionSteadyState(self, order=1):
        """This code computes the O(1/N) expansion of the mean field approximation
        (the term "V" is the "V" of Theorem~1 of https://hal.inria.fr/hal-01622054/document.

        Returns:
        if order = 0
            pi - fixed point of the mean field approximation
        if order = 1
            pi - fixed point of the mean field approximation
            V - fixed point of the refinement term
            (V, W) - fixed point values of V, W
        """
        pi = self.ode(time=10000)[1][-1]

        if order == 0:
            return pi
        if order >= 1:  # We need 2 derivatives and Q to get the O(1/N)-term
            Fp, Fpp, Q = self.defineDriftDerivativeQ(pi)

            # reshaping tensors
            N = self._model_dimension_N
            S = self._model_dimension_S
            pi = pi.flatten()
            Q = Q.reshape((N * S, N * S))
            Fp = Fp.reshape((N * S, N * S))
            Fpp = Fpp.reshape((N * S, N * S, N * S))

            # reduce dimension for steady state calculation
            Fp, Fpp, Q, P, Pinv, rank = self.reduceDimensionFpFppQ(Fp, Fpp, Q)
            if order == 1:
                pi, V, (V, W) = self.computePiV(pi, Fp, Fpp, Q)
                V, W = self.expandDimensionVW(V, W, Pinv)

                # reshaping into original representation
                pi = pi.reshape((N, S))
                V = V.reshape((N, S))
                W = W.reshape((N, S, N, S))
                return pi, V, (V, W)

    def reduceDimensionFpFppQ(self, Fp, Fpp, Q):
        P, P_inv, rank = self.dimensionReduction(Fp)
        Fp = (P @ Fp @ P_inv)[0:rank, 0:rank]
        Fpp = np.tensordot(np.tensordot(np.tensordot(P, Fpp, 1), P_inv, 1), P_inv,
                           axes=[[1], [0]])[0:rank, 0:rank, 0:rank]
        Q = (P @ Q @ P.transpose())[0:rank, 0:rank]
        return Fp, Fpp, Q, P, P_inv, rank

    def dimensionReduction(self, A):
        n = len(A)

        # M = np.array([l for l in self._list_of_transitions])
        M = self.get_transition_matrix()
        rank_of_transitions = np.linalg.matrix_rank(M)
        # If rank_of_transitions < n, this means that the stochastic process
        # evolves on a linear subspace of R^n.
        eigenvaluesOfJacobian = scipy.linalg.eig(A, left=False, right=False)
        rank_of_jacobian = np.linalg.matrix_rank(A)
        if sum(np.real(eigenvaluesOfJacobian) < 1e-8) < rank_of_transitions:
            # This means that there are less than "rank_of_transitions"
            # eigenvalues with <0 real part
            print("The Jacobian seems to be not Hurwitz")
        if rank_of_jacobian == n:
            return np.eye(n), np.eye(n), n
        C = np.zeros((n, n))
        n = len(A)
        d = 0
        rank_of_previous_submatrix = 0
        for i in range(n):
            rank_of_next_submatrix = np.linalg.matrix_rank(A[0:i + 1, 0:i + 1])
            if rank_of_next_submatrix > rank_of_previous_submatrix:
                C[d, i] = 1
                d += 1
            rank_of_previous_submatrix = rank_of_next_submatrix
        U, s, V = scipy.linalg.svd(A)
        C[rank_of_jacobian:, :] = U.transpose()[rank_of_jacobian:, :]
        return C, scipy.linalg.inv(C), rank_of_jacobian

    def get_transition_matrix(self):
        """Returns a Matrix with generated by appending all transitions of the system."""

        M = []

        trans_indices_A = np.nonzero(self.A)
        nr_unit_trans_A = trans_indices_A[0].shape[0]

        trans_indices_B = np.nonzero(self.B)
        nr_pair_trans_B = trans_indices_B[0].shape[0]
        # arrivals
        for i in range(nr_unit_trans_A):
            # transition indices
            k, s, s_prime = trans_indices_A[0][i], trans_indices_A[1][i], trans_indices_A[2][i]
            # add transition to M
            M.append((- self.e(k, s) + self.e(k, s_prime)).flatten())

        # pairwise transitions
        for i in range(nr_pair_trans_B):
            # transition indices
            k, k1 = trans_indices_B[0][i], trans_indices_B[1][i]
            s, s1 = trans_indices_B[2][i], trans_indices_B[3][i]
            s_prime, s1_prime = trans_indices_B[4][i], trans_indices_B[5][i]
            # add transition to M
            M.append((- self.e(k, s) + self.e(k, s_prime) - self.e(k1, s1) + self.e(k1, s1_prime)).flatten())

        return np.array(M)

    def computePiV(self, pi, Fp, Fpp, Q):
        """Returns the constants V and W (1/N-term for the steady-state)

        This function assumes that Fp is invertible.
        """
        from scipy.linalg import solve_continuous_lyapunov, inv

        # W = computeW(Fp, Q)
        W = solve_continuous_lyapunov(Fp, -Q)
        # V = computeV(Fp, Fpp, W)
        V = -np.tensordot(inv(Fp),
                          np.tensordot(Fpp, W / 2, 2),
                          1)
        return pi, V, (V, W)

    def expandDimensionVW(self, V, W, P_inv):
        rank = len(V)
        return P_inv[:, 0:rank] @ V, P_inv[:, 0:rank] @ W @ P_inv.transpose()[0:rank, :]

    def ode(self, time, number_of_steps=1000):
        """Simulates the ODE (mean-field approximation) up to a given time.

        Returns:
            T - numpy array specifying the times of the ode solution
            X - numpy array with values ( has shape T x N x S )
        """
        if self._x0 is None:
            raise InitialConditionNotDefined

        T = np.linspace(0, time, number_of_steps)

        def vector_valued_drift(x):
            # reshape input vector
            x = x.reshape((self._model_dimension_N, self._model_dimension_S))
            # return drift at x
            return self.defineDrift(x).flatten()

        X = integrate.odeint(lambda x, t: vector_valued_drift(x), self._x0.flatten(), T)

        X = np.array([X[i, :].reshape((self._model_dimension_N, self._model_dimension_S))
                      for i in range(X.shape[0])])

        return T, X

    def fixed_point(self):
        """Computes the fixed of the ODE (if this ODE has a fixed point starting from x0)."""
        if self._x0 is None:
            print(
                'No initial condition given. We assume that the initial condition is "x0=[[1,0,...],[1,0,...],...]"')
            self._x0 = np.zeros(self._model_dimension)
            self._x0[:, 0] = 1
        return super().fixed_point()

    def e(self, k, s):
        """Defines a unit vector with having the size of the model."""
        e_ = np.zeros(shape=(self._model_dimension_N, self._model_dimension_S))
        e_[k, s] = 1
        return e_
