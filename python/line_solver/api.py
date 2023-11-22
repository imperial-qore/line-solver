import jpype
import jpype.imports
from jpype import JInt

from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def cache_mva(gamma, m):
    ret = jpype.JPackage('jline').api.CACHE.cache_mva(jlineMatrixFromArray(gamma), jlineMatrixFromArray(m))
    pi = jlineMatrixFromArray(ret.pi)
    pi0 = jlineMatrixFromArray(ret.pi0)
    pij = jlineMatrixFromArray(ret.pij)
    x = jlineMatrixFromArray(ret.x)
    u = jlineMatrixFromArray(ret.u)
    E = jlineMatrixFromArray(ret.E)
    return pi, pi0, pij, x, u, E


def cache_prob_asy(gamma, m):
    return jpype.JPackage('jline').api.CACHE.cache_prob_asy(jlineMatrixFromArray(gamma), jlineMatrixFromArray(m))


def ctmc_uniformization(pi0, Q, t):
    return jlineMatrixToArray(
        jpype.JPackage('jline').api.CTMC.ctmc_uniformization(jlineMatrixFromArray(pi0), jlineMatrixFromArray(Q), t))


def ctmc_timereverse(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_timereverse(jlineMatrixFromArray(matrix)))


def ctmc_makeinfgen(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_makeinfgen(jlineMatrixFromArray(matrix)))


def ctmc_solve(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_solve(jlineMatrixFromArray(matrix)))


def ctmc_timereverse(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_timereverse(jlineMatrixFromArray(matrix)))


def dtmc_solve(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.DTMC.dtmc_solve(jlineMatrixFromArray(matrix)))


def dtmc_stochcomp(matrix, indexes):
    ind = jpype.java.util.ArrayList()
    for i in range(len(indexes)):
        ind.add(JInt(indexes[i]))
    return jlineMatrixToArray(jpype.JPackage('jline').api.DTMC.dtmc_stochcomp(jlineMatrixFromArray(matrix), ind))


def dtmc_timereverse(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.DTMC.dtmc_timereverse(jlineMatrixFromArray(matrix)))
