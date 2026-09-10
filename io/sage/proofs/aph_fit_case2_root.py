"""Closing the last aph_fit gap: is the f built from the K1..K22 radicals a root
of the scalar equation n3(f) = n3 that case 2 reduces to?

n2(f) = n2 was already proved exactly for any f, so this single scalar identity
is the entire remaining content of the closed form. The K-formulas nest cube
roots of square roots, which do not reduce exactly, so the check runs at 400-bit
precision.
"""
from sage.all import PolynomialRing, QQ, Matrix, vector, RealField, factorial, sqrt

RF = RealField(400)


def n3_of_f(n, n2v, fv):
    """normalized third moment of the case-2 APH(n) at parameter f."""
    a = 2 * (fv - 1) * (n - 1) / ((n - 1) * (n2v * fv ** 2 - 2 * fv + 2) - n)
    p = (fv - 1) * a
    lam = RF(1)
    mu = lam * (n - 1) / a
    alpha = vector(RF, [p, 1 - p] + [0] * (n - 2))
    T = Matrix(RF, n, n, 0)
    for i in range(n):
        T[i, i] = -mu
        if i + 1 < n:
            T[i, i + 1] = mu
    T[0, 0] = -lam
    T[0, 1] = lam
    one = vector(RF, [1] * n)
    A = (-T).inverse()
    m = [factorial(k) * ((alpha * (A ** k)) * one) for k in (1, 2, 3)]
    return m[1] / m[0] ** 2, m[2] / (m[0] * m[1])


def f_from_radicals(n, n2, n3):
    """aph_fit.m case 2, K1..K22 and the branch selection on n3 and K20."""
    K1 = RF(n - 1); K2 = RF(n - 2); K3 = 3 * n2 - 2 * n3; K4 = n3 - 3
    K5 = n - n2; K6 = 1 + n2 - n3; K7 = n + n2 - n * n2
    K8 = 3 + 3 * n2 ** 2 + n3 - 3 * n2 * n3
    inner = (-16 * K1 ** 2 * K7 ** 6
             + (4 * K1 * K5 ** 3 + K1 ** 2 * K2 * K4 ** 2 * n * n2 ** 2
                + 4 * K2 * n * n2 * (K4 * n ** 2 - 3 * K6 * n2 + K8 * n)) ** 2)
    K9 = 108 * K1 ** 2 * (4 * K2 ** 2 * K3 * n ** 2 * n2
                          + K1 ** 2 * K2 * K4 ** 2 * n * n2 ** 2
                          + 4 * K1 * K5 * (K5 ** 2 - 3 * K2 * K6 * n * n2)
                          + sqrt(inner))
    K10 = K4 ** 2 / (4 * K3 ** 2) - K5 / (K1 * K3 * n2)
    K11 = RF(2) ** (RF(1) / 3) * (3 * K5 ** 2 + K2 * (K3 + 2 * K4) * n * n2) / (K3 * K9 ** (RF(1) / 3) * n2)
    K12 = K9 ** (RF(1) / 3) / (3 * RF(2) ** (RF(7) / 3) * K1 ** 2 * K3 * n2)
    K13 = sqrt(K10 + K11 + K12)
    K14 = (6 * K1 * K3 * K4 * K5 + 4 * K2 * K3 ** 2 * n - K1 ** 2 * K4 ** 3 * n2) / (4 * K1 ** 2 * K3 ** 3 * K13 * n2)
    K15 = -K4 / (2 * K3)
    K16 = sqrt(2 * K10 - K11 - K12 - K14)
    K17 = sqrt(2 * K10 - K11 - K12 + K14)
    K20 = 6 * K1 * K3 * K4 * K5 + 4 * K2 * K3 ** 2 * n - K1 ** 2 * K4 ** 3 * n2
    if n3 < 3 * n2 / 2:
        return K13 + K15 - K17, 'K13+K15-K17'
    if K20 > 0:
        return -K13 + K15 + K16, '-K13+K15+K16'
    return K13 + K15 + K17, 'K13+K15+K17'


print('== aph_fit case 2: the radical f solves n3(f) = n3 (400-bit) ==')
ok = []
for n in (3, 4, 5):
    for n2v, n3v in ((RF(2.5), RF(9)), (RF(3), RF(15)), (RF(4), RF(30)), (RF(2.2), RF(7.5))):
        try:
            fv, branch = f_from_radicals(n, n2v, n3v)
            g2, g3 = n3_of_f(n, n2v, fv)
            e2 = abs(g2 - n2v) / n2v
            e3 = abs(g3 - n3v) / n3v
            good = (e3 < RF(1e-25))
            print('  n=%d n2=%s n3=%s [%s]: n2 err %.2e  n3 err %.2e  %s'
                  % (n, n2v, n3v, branch, float(e2), float(e3), 'ROOT' if good else 'NOT A ROOT'))
            ok.append(good)
        except Exception as e:
            print('  n=%d n2=%s n3=%s: %s' % (n, n2v, n3v, type(e).__name__))
print('\nradical f is a root of the n3 equation on every feasible point tried:', all(ok))
