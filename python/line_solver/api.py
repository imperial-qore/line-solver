import jpype
import jpype.imports

from line_solver import jlineMatrixToArray, jlineMatrixFromArray

def ctmc_uniformization(pi0, Q, t):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_uniformization(jlineMatrixFromArray(pi0), jlineMatrixFromArray(Q), t))
def ctmc_timereverse(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_timereverse(jlineMatrixFromArray(matrix)))
def ctmc_makeinfgen(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_makeinfgen(jlineMatrixFromArray(matrix)))
def ctmc_solve(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.CTMC.ctmc_solve(jlineMatrixFromArray(matrix)))
def dtmc_solve(matrix):
    return jlineMatrixToArray(jpype.JPackage('jline').api.DTMC.dtmc_solve(jlineMatrixFromArray(matrix)))