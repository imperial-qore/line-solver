"""Pre-compute the MMAP(2,K) marking inverse in closed form.

The 3K x 3K system is block diagonal: class c's fractions (q1c, q2c, q3c)
affect only (p_c, p_c F_c, p_c B_c). So the inverse is ONE 3x3 block, identical
for every class and independent of K. What follows prints that block, and the
resulting per-class formulas, for both canonical forms, ready to be hard-coded.
"""
from sage.all import PolynomialRing, QQ, Matrix, vector

for form in (1, 2):
    R = PolynomialRing(QQ, ['h1', 'h2', 'r1', 'r2', 'q1', 'q2', 'q3'], order='lex')
    F = R.fraction_field()
    h1, h2, r1, r2, q1, q2, q3 = [F(g) for g in R.gens()]
    D0 = Matrix(F, [[-1 / h1, r1 / h1], [0, -1 / h2]])
    if form == 1:
        D1 = Matrix(F, [[(1 - r1) / h1, 0], [(1 - r2) / h2, r2 / h2]])
        D1c = Matrix(F, [[D1[0, 0] * q1, 0], [D1[1, 0] * q2, D1[1, 1] * q3]])
    else:
        D1 = Matrix(F, [[0, (1 - r1) / h1], [(1 - r2) / h2, r2 / h2]])
        D1c = Matrix(F, [[0, D1[0, 1] * q1], [D1[1, 0] * q2, D1[1, 1] * q3]])

    one = vector(F, [1, 1])
    A = (-D0).inverse()
    P = A * D1
    T = (P.transpose() - Matrix.identity(F, 2)); T[1, :] = Matrix(F, 1, 2, [1, 1])
    pie = T.solve_right(vector(F, [0, 1]))

    pc = (pie * (A * D1c)) * one
    pB = (pie * (A * A * D1c)) * one          # = p_c * B_c
    pF = (pie * (A * D1c * A)) * one          # = p_c * F_c
    y = [pc, pF, pB]
    qs = [q1, q2, q3]
    M = Matrix(F, [[F(F(y[i]).derivative(R.gen(4 + j))) for j in range(3)] for i in range(3)])
    c0 = vector(F, [F(y[i]) - sum(M[i, j] * qs[j] for j in range(3)) for i in range(3)])

    print('=' * 70)
    print('canonical form %d' % form)
    print('  affine offset c0 =', [F(x) for x in c0])
    print('  det M =', F(M.determinant()).factor())
    Mi = M.inverse()
    print('  inverse block M^-1 (rows give q1, q2, q3 from [p_c, p_c F_c, p_c B_c]):')
    for i in range(3):
        for j in range(3):
            e = F(Mi[i, j])
            print('    Minv[%d][%d] = %s' % (i, j, e.factor() if e != 0 else '0'))
    # the resulting per-class closed form, in the shape the m3a code uses
    print('  => q_j = sum_i Minv[j][i] * y_i   with y = (p_c, p_c F_c, p_c B_c)')
    print('     equivalently q_j = p_c * (Minv[j][0] + Minv[j][1] F_c + Minv[j][2] B_c)')
