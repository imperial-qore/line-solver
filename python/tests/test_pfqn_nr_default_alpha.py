"""The Norlund-Rice rate matrix defaults to the load-independent mu(i,n) = 1.

`alpha` is the (M x Ntot) load-dependent rate matrix `mu`: a single-server queue
is mu(i,n) = 1 and a delay is mu(i,n) = n, which is why pfqn_ncld appends 1:Ntot
as the think-time ROW and leaves every other row alone. pfqn_nrl/pfqn_nrp used to
default the WHOLE matrix to 1:Ntot, making every station an infinite server, so
pfqn_nc -- the load-independent dispatcher, and the only caller that omits alpha
-- returned the constant of a different model without saying so. pfqn_nre always
defaulted to ones, so nre alone was right and the two wrong methods looked
merely inaccurate rather than misdirected.
"""

import numpy as np

from line_solver.api.pfqn import pfqn_nre, pfqn_nrl, pfqn_nrp
from line_solver.api.pfqn.nc import pfqn_ca, pfqn_nc
from line_solver.api.pfqn.ncld import pfqn_ncld

# three single-server stations, two classes; the model of test_pfqn_nre
L = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
N = np.array([2.0, 3.0])
NTOT = 5
ONES = np.ones((3, NTOT))
ALL_DELAY = np.tile(np.arange(1.0, NTOT + 1), (3, 1))

LG_EXACT_NOZ = 4.139104712019314      # pfqn_ca, exact for this model class
LG_EXACT_Z = 4.783752287503922        # with Z = [1, 0.5]; matches test_pfqn_nre


def test_omitted_alpha_is_the_load_independent_rate_matrix():
    Z = np.zeros(2)
    for fn in (pfqn_nrl, pfqn_nrp):
        assert fn(L, N, Z) == fn(L, N, Z, alpha=ONES)
    assert abs(pfqn_nre(L, N, Z) - pfqn_nre(L, N, Z, alpha=ONES)) < 1e-9


def test_pfqn_nc_matches_the_load_dependent_dispatcher_at_mu_one():
    # what SolverNC computes for these three names on a load-independent model
    Z = np.zeros(2)
    for method in ('nrl', 'nrp', 'nre'):
        _, lG = pfqn_nc(L, N, Z, method=method)
        assert abs(lG - pfqn_ncld(L, N, Z, ONES, {'method': method}).lG) < 1e-9


def test_the_all_delay_matrix_is_a_different_model():
    # the guard on the regression: mu(i,n) = n is a valid rate matrix, it just
    # describes a network of infinite servers, whose constant is nowhere near
    Z = np.zeros(2)
    assert abs(pfqn_ca(L, N, Z)[1] - LG_EXACT_NOZ) < 1e-12
    for method in ('nrl', 'nrp', 'nre'):
        li = pfqn_ncld(L, N, Z, ONES, {'method': method}).lG
        delay = pfqn_ncld(L, N, Z, ALL_DELAY, {'method': method}).lG
        assert abs(li - LG_EXACT_NOZ) < 0.5
        assert delay < LG_EXACT_NOZ - 2.0


def test_think_time_row_rides_on_top_of_the_ones_default():
    # Z > 0 appends 1:Ntot as ONE more row; the queueing rows stay at one
    Z = np.array([1.0, 0.5])
    assert abs(pfqn_ca(L, N, Z)[1] - LG_EXACT_Z) < 1e-12
    for method in ('nrl', 'nrp', 'nre'):
        _, lG = pfqn_nc(L, N, Z, method=method)
        assert abs(lG - pfqn_ncld(L, N, Z, ONES, {'method': method}).lG) < 1e-9
        assert abs(lG - LG_EXACT_Z) < 0.5
