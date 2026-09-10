"""SageMath proof of the m3a marking-coefficient identities.

maph2m_fit_multiclass and mamap2m_fit_fb_multiclass both claim that the class
marking fractions q(j,c) are AFFINE in the class characteristics:

    MAPH:  q(j,c) = B(c) * q_b(j,c) + q_0(j,c)
    MAMAP: q(j,c) = F(c) * q_f(j,c) + B(c) * q_b(j,c) + q_0(j,c)

with coefficients that depend only on the underlying (h1, h2, r1, r2) and the
class probability p(c). Those affine relations are what turn the fit into a
quadratic program; if a coefficient is wrong the fit stays self-consistent and
a round-trip test cannot see it.

Here the marking fractions q are the FREE variables: p(c), F(c) and B(c) are
computed from the MMAP with the m3a definitions (mmap_pc.m,
mmap_forward_moment.m, mmap_backward_moment.m), substituted into the claimed
right-hand side, and the result must reduce to q identically.
"""
from sage.all import PolynomialRing, QQ, Matrix, vector

R = PolynomialRing(QQ, ['h1', 'h2', 'r1', 'r2', 'q1', 'q2', 'q3'], order='lex')
F = R.fraction_field()
h1, h2, r1, r2, q1, q2, q3 = [F(g) for g in R.gens()]


def pie_of(D0, D1):
    """map_pie: stationary vector of the embedded DTMC P = (-D0)^-1 D1."""
    n = D0.nrows()
    P = (-D0).inverse() * D1
    M = (P.transpose() - Matrix.identity(F, n)).augment(Matrix(F, n, 0))
    M = M.stack(Matrix(F, 1, n, [1] * n))
    # solve pi (P - I) = 0 with sum(pi) = 1 by replacing the last row
    A = (P.transpose() - Matrix.identity(F, n))
    A[n - 1, :] = Matrix(F, 1, n, [1] * n)
    rhs = vector(F, [0] * (n - 1) + [1])
    return A.solve_right(rhs)


def class_stats(D0, D1, D1c):
    """(p_c, F_c, B_c) with the m3a definitions."""
    n = D0.nrows()
    one = vector(F, [1] * n)
    pie = pie_of(D0, D1)
    M = (-D0).inverse()
    pc = (pie * (M * D1c)) * one
    # backward: 1!/pc * sum(pie * M^2 * D1c)
    Bc = ((pie * (M * M * D1c)) * one) / pc
    # forward: 1!/pc * sum(pie * (M D1c) * M)
    Fc = ((pie * (M * D1c * M)) * one) / pc
    return pc, Fc, Bc


def report(name, expr):
    ok = (F(expr) == 0)
    print('  %-56s %s' % (name, 'IDENTITY' if ok else 'NOT ZERO: %s' % F(expr)))
    return ok


results = []

print('== MAPH(2,m): q(j,c) = B(c) q_b(j,c) + q_0(j,c) ==')
# canonical acyclic APH(2): amap2_assemble(h1, h2, r1, 0, form 1)
D0 = Matrix(F, [[-1 / h1, r1 / h1], [0, -1 / h2]])
D1 = Matrix(F, [[(1 - r1) / h1, 0], [1 / h2, 0]])
D1c = Matrix(F, [[D1[0, 0] * q1, 0], [D1[1, 0] * q2, 0]])
pc, Fc, Bc = class_stats(D0, D1, D1c)
qb1 = pc * (1 / (h2 * (r1 - 1)))
q01 = pc * (-(h1 + h2) / (h2 * (r1 - 1)))
qb2 = pc * (1 / (h2 * r1))
q02 = pc * (-h1 / (h2 * r1))
results.append(report('q1 - (B q_b(1) + q_0(1))', q1 - (Bc * qb1 + q01)))
results.append(report('q2 - (B q_b(2) + q_0(2))', q2 - (Bc * qb2 + q02)))

print('== MAMAP(2,m) F+B, first canonical form ==')
# amap2_assemble(h1, h2, r1, r2, form 1); marking mask [q1 0; q2 q3]
D0 = Matrix(F, [[-1 / h1, r1 / h1], [0, -1 / h2]])
D1 = Matrix(F, [[(1 - r1) / h1, 0], [(1 - r2) / h2, r2 / h2]])
D1c = Matrix(F, [[D1[0, 0] * q1, 0], [D1[1, 0] * q2, D1[1, 1] * q3]])
pc, Fc, Bc = class_stats(D0, D1, D1c)
qf = [F(0),
      -(pc * (r1 * r2 - r2 + 1)) / (r1 * (h1 + h2 * (r1 - 1)) * (r2 - 1)),
      -(pc * (r1 * r2 - r2 + 1)) / (r1 * r2 * (h1 - h2 + h2 * r1))]
qb = [-(pc * (r1 * r2 - r2 + 1)) / ((h2 - h1 * r2) * (r1 - 1) * (r2 - 1)),
      -(pc * (r1 * r2 - r2 + 1)) / (r1 * (h2 - h1 * r2) * (r2 - 1)),
      F(0)]
q0 = [(pc * (h1 + h2 - h1 * r2) * (r1 * r2 - r2 + 1)) / ((h2 - h1 * r2) * (r1 - 1) * (r2 - 1)),
      ((pc * (r1 * r2 - r2 + 1)) / ((r1 - 1) * (r2 - 1))
       + (h1 * pc * (r1 * r2 - r2 + 1)) / (r1 * (h2 - h1 * r2) * (r2 - 1))
       - (h1 * pc * (r1 * r2 - r2 + 1)) / (r1 * (h1 + h2 * (r1 - 1)) * (r1 - 1) * (r2 - 1))),
      (pc * (h1 + h2 * r1) * (r1 * r2 - r2 + 1)) / (r1 * r2 * (h1 - h2 + h2 * r1))]
for j, qj in enumerate([q1, q2, q3]):
    results.append(report('form 1: q%d - (F q_f + B q_b + q_0)' % (j + 1),
                          qj - (Fc * qf[j] + Bc * qb[j] + q0[j])))

print('== MAMAP(2,m) F+B, second canonical form ==')
# amap2_assemble(h1, h2, r1, r2, form 2); marking mask [0 q1; q2 q3]
D0 = Matrix(F, [[-1 / h1, r1 / h1], [0, -1 / h2]])
D1 = Matrix(F, [[0, (1 - r1) / h1], [(1 - r2) / h2, r2 / h2]])
D1c = Matrix(F, [[0, D1[0, 1] * q1], [D1[1, 0] * q2, D1[1, 1] * q3]])
pc, Fc, Bc = class_stats(D0, D1, D1c)
qf = [F(0),
      (pc * (r1 + r2 - r1 * r2 - 2)) / ((r2 - 1) * (h1 - h2 + h2 * r1)),
      (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 + h2 * (r1 - 1)))]
qb = [-(pc * (r1 + r2 - r1 * r2 - 2)) / ((r1 - 1) * (r2 - 1) * (h1 - h2 - h1 * r1 + h1 * r1 * r2)),
      F(0),
      (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2))]
q0 = [(pc * (h2 + h1 * r1 - h1 * r1 * r2) * (r1 + r2 - r1 * r2 - 2)) / ((r1 - 1) * (r2 - 1) * (h1 - h2 - h1 * r1 + h1 * r1 * r2)),
      -(h2 * pc * (r1 + r2 - r1 * r2 - 2)) / ((r2 - 1) * (h1 - h2 + h2 * r1)),
      ((h1 * pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 + h2 * (r1 - 1)) * (r1 - 1))
       - (h1 * pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2))
       - (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (r1 - 1)))]
for j, qj in enumerate([q1, q2, q3]):
    results.append(report('form 2: q%d - (F q_f + B q_b + q_0)' % (j + 1),
                          qj - (Fc * qf[j] + Bc * qb[j] + q0[j])))

print('== MAMAP(2,m) F+B, non-canonical phase-type branch (r1 = 1) ==')
# same first canonical form at r1 = 1, where only the forward moments are fitted
Rn = PolynomialRing(QQ, ['h1', 'h2', 'r2', 'q2', 'q3'], order='lex')
Fn = Rn.fraction_field()
n_h1, n_h2, n_r2, n_q2, n_q3 = [Fn(g) for g in Rn.gens()]
D0n = Matrix(Fn, [[-1 / n_h1, 1 / n_h1], [0, -1 / n_h2]])
D1n = Matrix(Fn, [[0, 0], [(1 - n_r2) / n_h2, n_r2 / n_h2]])
D1cn = Matrix(Fn, [[0, 0], [D1n[1, 0] * n_q2, D1n[1, 1] * n_q3]])
onen = vector(Fn, [1, 1])
Pn = (-D0n).inverse() * D1n
An = (Pn.transpose() - Matrix.identity(Fn, 2))
An[1, :] = Matrix(Fn, 1, 2, [1, 1])
pien = An.solve_right(vector(Fn, [0, 1]))
Mn = (-D0n).inverse()
pcn = (pien * (Mn * D1cn)) * onen
Fcn = ((pien * (Mn * D1cn * Mn)) * onen) / pcn
r1n = Fn(1)
qf_nc = [pcn * (-1 / ((n_h1 + n_h2 * (r1n - 1)) * (n_r2 - 1) * (r1n + n_r2 - r1n * n_r2))),
         pcn * (-1 / (n_r2 * (n_h1 + n_h2 * (r1n - 1)) * (r1n + n_r2 - r1n * n_r2)))]
q0_nc = [pcn * (n_h2 / ((n_r2 - 1) * (r1n + n_r2 - r1n * n_r2) * (n_h1 - n_h2 + n_h2 * r1n))),
         pcn * ((n_h1 + n_h2 * r1n) / (n_r2 * (r1n + n_r2 - r1n * n_r2) * (n_h1 - n_h2 + n_h2 * r1n)))]
for j, qj in enumerate([n_q2, n_q3]):
    val = qj - (Fcn * qf_nc[j] + q0_nc[j])
    ok = (Fn(val) == 0)
    print('  %-56s %s' % ('r1=1: q%d - (F q_f + q_0)' % (j + 2),
                          'IDENTITY' if ok else 'NOT ZERO: %s' % Fn(val)))
    results.append(ok)

print('\nall marking-coefficient identities hold:', all(results))
