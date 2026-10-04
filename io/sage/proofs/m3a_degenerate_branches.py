"""Degenerate-branch marking coefficients of mamap2m_fit_fb_multiclass (form 2,
r2 = 0), where only the forward OR the backward moments can be fitted."""
from sage.all import PolynomialRing, QQ, Matrix, vector

R = PolynomialRing(QQ, ['h1', 'h2', 'r1', 'q1', 'q2'], order='lex')
F = R.fraction_field()
h1, h2, r1, q1, q2 = [F(g) for g in R.gens()]

# second canonical form with r2 = 0
D0 = Matrix(F, [[-1 / h1, r1 / h1], [0, -1 / h2]])
D1 = Matrix(F, [[0, (1 - r1) / h1], [1 / h2, 0]])
D1c = Matrix(F, [[0, D1[0, 1] * q1], [D1[1, 0] * q2, 0]])

one = vector(F, [1, 1])
P = (-D0).inverse() * D1
T = (P.transpose() - Matrix.identity(F, 2))
T[1, :] = Matrix(F, 1, 2, [1, 1])
pie = T.solve_right(vector(F, [0, 1]))
M = (-D0).inverse()
pc = (pie * (M * D1c)) * one
Bc = ((pie * (M * M * D1c)) * one) / pc
Fc = ((pie * (M * D1c * M)) * one) / pc

res = []
print('== form 2, r2 = 0: forward variant ==')
qf = [pc * (-(r1 - 2) / ((h1 + h2 * (r1 - 1)) * (r1 - 1))),
      pc * (-(r1 - 2) / (h1 + h2 * (r1 - 1)))]
q0 = [pc * (1 - (h1 + h2) / ((r1 - 1) * (h1 - h2 + h2 * r1))),
      pc * ((h2 * (r1 - 2)) / (h1 - h2 + h2 * r1))]
for j, qj in enumerate([q1, q2]):
    val = qj - (Fc * qf[j] + q0[j])
    ok = (F(val) == 0)
    print('  q%d - (F q_f + q_0): %s' % (j + 1, 'IDENTITY' if ok else 'NOT ZERO: %s' % F(val)))
    res.append(ok)

print('== form 2, r2 = 0: backward variant ==')
qb = [pc * (-(r1 - 2) / ((h2 + h1 * (r1 - 1)) * (r1 - 1))),
      pc * (-(r1 - 2) / (h2 + h1 * (r1 - 1)))]
q0b = [pc * (1 - (h1 + h2) / ((r1 - 1) * (h2 - h1 + h1 * r1))),
       pc * ((h1 * (r1 - 2)) / (h2 - h1 + h1 * r1))]
for j, qj in enumerate([q1, q2]):
    val = qj - (Bc * qb[j] + q0b[j])
    ok = (F(val) == 0)
    print('  q%d - (B q_b + q_0): %s' % (j + 1, 'IDENTITY' if ok else 'NOT ZERO: %s' % F(val)))
    res.append(ok)

print('\ndegenerate-branch coefficients hold:', all(res))
