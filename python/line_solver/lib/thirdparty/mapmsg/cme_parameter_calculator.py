"""
CME Parameter Calculator.

Port of CMEParameterCalculator.m from the MAPMsG MATLAB library.

Computes Concentrated Matrix Exponential (CME) parameters for approximating
deterministic distributions, used in first passage time calculations.

Reference:
    O. Gursoy, K. A. Mehr, N. Akar, "The MAP/M/s + G Call Center Model with
    Generally Distributed Patience Times."
"""

import numpy as np
from numpy.linalg import inv


# ---------------------------------------------------------------------------
# Pre-computed CME data tables
# ---------------------------------------------------------------------------
# b vectors (initial probability weights) and off-diagonal step values for
# the block-structured B matrices, indexed by (n-1)/2.

_b_data = {
    12: [
        201.29305527432533, 121.95756811242167, -193.096102096128,
        -122.68122499055721, -81.69176197462747, -63.22266724705794,
        65.57763017385028, 28.89925358153127, 46.53109839331482,
        30.299758102542366, -8.738691659155492, -0.001179508655788386,
        -17.045775482071868, -8.155939777086843, -2.2790187663215016,
        -1.888556205662306, 3.2583909875568158, 1.0603540479937295,
        0.9809562591378991, 0.35863084865038025, -0.2717496542617506,
        -0.05183626531614445, -0.08736194837838565, -0.010901997305688587,
        0.0060717912617999316,
    ],
    25: [
        1111.3459597426531, 1013.5935750153698, -906.7724040373613,
        -482.6223341096786, -888.5363985989719, -736.5413252484491,
        94.61104028208489, -172.52385336054593, 528.8635023684683,
        296.5461156633309, 303.6649723656398, 315.8841492300176,
        -87.12864980829919, 62.91229008242268, -246.45535524687156,
        -140.35453625111003, -138.7964149773478, -147.49931732158603,
        38.03392139557178, -33.63626666074197, 111.30029777975795,
        57.964748933608966, 65.32582650776743, 63.28943685171059,
        -10.885363786007277, 17.348580959997935, -42.5907221256393,
        -18.737356341772372, -25.822464844580605, -21.261270824071012,
        1.6032761370814057, -6.1778960545226695, 12.161620791711542,
        4.450757575057735, 6.997585731699266, 4.74259587045975,
        -0.22615038425938033, 1.1015824998685881, -2.2774448501728464,
        -0.746303694292432, -0.9833727638335755, -0.5268272080403923,
        0.11401418980312379, -0.03592806724472174, 0.20000686706058754,
        0.05436432930524405, 0.030428482032967934, 0.009614779828734247,
        -0.009966881166604605, -0.0010028568589930972, -0.001338128885418459,
    ],
    50: [
        6402.224420062449, 7110.574657198787, -3977.1390106717413,
        -681.0341376071109, -6443.20793478582, -4632.445983450395,
        -2054.5155786722044, -3426.3862164706793, 2143.2834833155753,
        -252.35139857632677, 3290.2039217147785, 2035.975051654224,
        1846.5707481446416, 2296.3508259921828, -367.0100363905309,
        1002.3464307521068, -1720.2640260608136, -589.496943039846,
        -1615.6215147707965, -1405.2895157933654, -511.0465796972745,
        -1134.0108634626024, 648.3863858591234, -211.9554858273741,
        1125.0849135567719, 629.6650276520397, 783.9901425295026,
        883.7525103544706, 29.269856166402647, 527.0602188558557,
        -572.1473668064381, -77.67095262832146, -679.7883809392939,
        -496.22310940228425, -339.548295220095, -510.44234446121413,
        133.44392984665106, -205.24488303788644, 413.8643326292256,
        154.60571752605972, 372.57742661237194, 332.7214842357738,
        111.96678847102146, 262.9680160934082, -153.18810823938847,
        50.09247991795238, -257.8543309185617, -138.1432081056258,
        -178.38137989913918, -192.38277159813407, -11.880882725626307,
        -115.47907294741796, 116.00874620646528, 9.001784685314545,
        137.77886253958, 91.41043303810943, 70.71603765377363,
        94.19873617702176, -17.848539226283776, 40.420700059328375,
        -67.53889998904101, -19.086402276341683, -60.91914139766393,
        -46.50854899969573, -21.059524667187496, -36.74067928204501,
        16.253860799060572, -9.526473069558264, 29.53976494149302,
        11.915634713636276, 20.23711699477028, 17.04129758946868,
        3.300835328863558, 9.865808377178183, -7.670992414179801,
        0.44405887365957425, -8.724897788201377, -4.308563015312876,
        -4.04224524473775, -3.8174639796122074, 0.48464312169329876,
        -1.2551802498761244, 2.047090957869685, 0.5059169610722442,
        1.3331444739463298, 0.7692227463010747, 0.21522772872447127,
        0.32170752742924347, -0.2656245704979173, -0.02326024000224338,
        -0.19852674948065543, -0.07852019677306908, -0.034559206207524305,
        -0.024773714706152774, 0.019932498138524764, 0.002502779081987538,
        0.009239211765793214, 0.0020779771735582107, 0.00013822743623946983,
        0.00007132362879156872, -0.00023416964833831966,
    ],
}

_step_data = {
    12: 0.66509,
    25: 0.557478,
    50: 0.478865,
}


class MESystem:
    """Minimal state-space representation (A, B, C) for a matrix exponential distribution.

    Mirrors MATLAB's ``ss(A, B, C, 0)`` structure used in the original code.

    Attributes
    ----------
    A : ndarray, shape (n, n)
        State transition matrix.
    B : ndarray, shape (n, 1)
        Input (column) vector, equal to ``-A @ ones``.
    C : ndarray, shape (1, n)
        Output (row) vector (initial probability vector).
    """

    def __init__(self, A, B, C):
        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.C = np.asarray(C, dtype=float)


def _build_B_matrix(n, step):
    """Construct the n x n block-structured B matrix.

    Structure: first diagonal element is -1 (1x1 block), then (n-1)/2
    blocks of size 2x2::

        [-1,      -k*step]
        [ k*step, -1     ]

    for k = 1, 2, ..., (n-1)/2.
    """
    B = np.zeros((n, n))
    B[0, 0] = -1.0
    num_blocks = (n - 1) // 2
    for k in range(1, num_blocks + 1):
        idx = 1 + 2 * (k - 1)
        val = k * step
        B[idx, idx] = -1.0
        B[idx, idx + 1] = -val
        B[idx + 1, idx] = val
        B[idx + 1, idx + 1] = -1.0
    return B


def cme_parameter_calculator(n, time_delay):
    """Compute CME (Concentrated Matrix Exponential) parameters.

    Parameters
    ----------
    n : int
        Order of the CME approximation.  Must be 25, 51, or 101.
    time_delay : float
        The time horizon (tau) used to scale the CME distribution.

    Returns
    -------
    me_system : MESystem
        State-space representation ``(A, B, C)`` of the scaled CME distribution.
    cv_reciprocal : float
        Reciprocal of the squared coefficient of variation.

    Notes
    -----
    The MATLAB code maps ``OrderofPHCME`` values 1/2/3 to ``n`` = 25/51/101,
    then computes ``ord = (n - 1) / 2`` to index into the pre-computed tables
    (keyed by 12, 25, 50).
    """
    ord_val = (n - 1) // 2
    if ord_val not in _b_data:
        raise ValueError(
            f"Unsupported CME order n={n} (ord={ord_val}). "
            f"Supported orders: n=25 (ord=12), n=51 (ord=25), n=101 (ord=50)."
        )

    b8 = np.array(_b_data[ord_val])
    step = _step_data[ord_val]
    B8 = _build_B_matrix(n, step)

    # MATLAB: MESystem = ss(B8, B8*ones(n,1), -b8, 0)  [intermediate, unused]
    # Then:   m1 = -b8 * inv(B8) * ones(n,1)
    e = np.ones((n, 1))
    B8bar = inv(B8)
    m1_scalar = float(-b8 @ B8bar @ e)

    # MATLAB: MESystem = ss(B8*m1/TimeDelay, -B8*m1/TimeDelay*ones(n,1), b8, 0)
    A = B8 * m1_scalar / time_delay
    B_vec = -A @ e  # = -(B8 * m1 / TimeDelay) * ones = -A * e
    C = b8.reshape(1, -1)

    me_system = MESystem(A, B_vec, C)

    # Compute coefficient of variation
    Abar = inv(me_system.A)
    m1_check = float(me_system.C @ Abar @ Abar @ me_system.B)
    m2_check = float(-2.0 * me_system.C @ np.linalg.matrix_power(Abar, 3) @ me_system.B)
    cv = (m2_check - m1_check ** 2) / (m1_check ** 2)
    cv_reciprocal = 1.0 / cv

    return me_system, cv_reciprocal
