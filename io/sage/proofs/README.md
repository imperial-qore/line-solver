# Symbolic validation of the closed-form fitters

Each script here proves, in exact arithmetic, that one of the closed-form
fitting formulas in `matlab/lib/kpctoolbox` or `matlab/lib/m3a` really does what
it claims. They are the reason those formulas can be trusted beyond a numerical
round trip: a round trip only shows that a fitter is self-consistent, while a
mis-derived coefficient can be self-consistent and still wrong.

## Running them

They need SageMath, which the repository already ships as a container:

```bash
docker run --rm -i imperialqore/line-sage-rest sage -python - < proofs/mmpp2_fit3.py
```

The script is piped on stdin rather than mounted: a confined Docker cannot bind
mount an arbitrary scratch directory (see `_kb/11-conventions-and-gotchas.md`).
Runtimes are seconds to a couple of minutes each.

## Technique

Radicals are the only obstacle to exact reduction, and they are handled the same
way everywhere: the radical is replaced by an indeterminate `RAD` and every
identity is reduced modulo `RAD^2 - DISC` in a lex ring with `RAD` largest.
Writing the reduced numerator as `A + RAD*B`, the identity holds on BOTH
branches of the square root iff `A = B = 0`. Where the closed form inverts a
map (AMAP(2), APH(2), MAP(2)), the free variables are the PARAMETERS: the
characteristics are computed from them with the standard MAP formulas, fed into
the closed-form inverse, and the parameters must come back.

## What each script covers

| script | proves |
|---|---|
| `mmpp2_fit3.py` | the MMPP(2) rates reproduce E1, E2, E3, gamma2, rho1; plus `g2 = tr(P)-1`, `SCV >= 1`, `rho1 = (g2/2)(1-1/SCV)`, `rho_k` geometric, `g2 = (IDC-SCV)/(IDC-1)` |
| `map2_fit.py` | `map2_fit` reproduces (e1, e2, e3, g2) on all three branches (hyper, correlated, hypo) |
| `aph2_amap2_hyperexp.py` | `aph2_fitall` and `amap2_fitall_gamma` (both canonical forms) invert correctly; `map_hyperexp` matches mean and SCV on both roots and its reachable SCV at fixed p is (2-p)/p; `(3/2)E2^2/E1` is a strict floor on E3 over MMPP(2) |
| `m3a_marking_coefficients.py` | the MAPH and MAMAP marking fractions really are affine in the forward/backward moments, with the coefficients the fitters use (11 identities) |
| `m3a_degenerate_branches.py` | the same for the degenerate form-2 branches, forward and backward variants |
| `mmpp2_counting_idc.py` | derives Var[N(t)] for the MMPP(2) from first principles and confirms the IDC(t) shape `mmpp2_fitc` inverts through Lambert W, the closed form of IDC(inf), and the `xbt1` term of the `mmpp2_fitc_approx` objective |
| `mmpp2_counting_third_moment.py` | the same derivation to third order: the third central moment of counts equals the `xm3t2` term of that objective |
| `aph_fit_orders.py` | `aph_fit` case 1 matches n2 and n3 exactly for orders 2 to 6, and case 2 matches n2 for ANY value of its free parameter f |
| `qlen_tail_moments.py` | the survival (tail) vertex: the edge and its inverse against the definition with FREE masses, univariate and joint; the single-class survival identity `P(n >= k) = (prod L^k) G(N-|k|)/G(N)` with free demands and think time; the multiclass complementary-network law; and that the naive survival form FAILS for two classes, which is why `pfqn_qlen_joint_moments` has two routes |
| `moment_cumulants.py` | the cumulant vertex: with the moments left as FREE indeterminates, the implemented recursions equal the coefficients of log/exp of the exponential generating function, and both Leonov-Shiryaev set-partition formulas; the factorial cumulants against the log of the probability generating function, with the Poisson case |
| `moment_multivariate.py` | the joint conversions: the twelve separable edges on product point masses (d = 2 and 3), the joint central pair with FREE means, the joint cumulants against the joint cumulant generating function (FREE joint moments) plus the multivariate Leonov-Shiryaev formula, the multinomial marking rule from the composition of the pgf, and the aggregation rule from the Vandermonde convolution of falling factorials |
| `moment_conversions.py` | the whole Heindl-van de Liefvoort house of moments (`matlab/src/api/moment`, `jline.api.moment`, `line_solver.api.moment`): the four triangles against Sage's own Stirling/Lah, the twelve conversions against the DEFINITION of the six moment families, the central pair on a symbolic two-atom law, and the commutativity of every path in the house |
| `mmap2k_marking_inverse.py` | the MMAP(2,K) marking system is block diagonal, so its inverse is ONE 3x3 block, identical per class and independent of K; emits the closed form used by `mmap2k_fit` in all three codebases |
| `mmap2k_degrees_of_freedom.py` | the marked acyclic canonical form loses no generality: a general MMAP(2,K) and the canonical family both realize 3K+1 independent characteristics (Jacobian rank, K = 1..6) |
| `mmapnk_identifiability.py` | at orders 2 to 5 the span of every characteristic linear in the class matrix has full rank z = nnz(D1), so the marking is identifiable at any order; the MMAP(2,K) argument is not special to order two |
| `mmap3k_marking_inverse.py` | the independent characteristic set of lowest total order per order n -- (1,0),(1,1),(2,0) at n=2 and (1,0),(1,1),(2,0),(3,0) at n=3 -- with the size of the symbolic inverse that would have to be inlined |
| `aph_fit_case2_root.py` | the f built from the K1..K22 radicals is a root of the single scalar equation case 2 reduces to, at 400-bit precision |

## Not covered

- The FEASIBILITY BOUNDS (`amap2_adjust_gamma` gamma interval, `aph2_adjust`,
  the order search in `aph_fit`) are extremal characterizations from the
  literature, not algebraic identities; proving them needs an optimization
  argument, not a reduction.
- `mamap22_fit_fs_multiclass` / `mamap22_fit_bs_multiclass`, which are not
  ported to native Python and are unaudited in the JAR.
- The branch SELECTION rules of `aph_fit` case 2 (which root to take): the
  400-bit check confirms the returned f is a root of the right equation.
