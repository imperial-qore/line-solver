
from sage.all import PolynomialRing, QQ, Matrix, vector, factorial

R = PolynomialRing(QQ, ['RAD', 'E1', 'SCV', 'E3', 'G2'], order='lex')
RAD, E1, SCV, E3, G2 = R.gens()
F = R.fraction_field()
RAD, E1, SCV, E3, G2 = [F(g) for g in R.gens()]
DISCR = R('E3**2-12*E1**3*SCV*E3+6*E1**3*G2*E3-6*G2*SCV*E1**3*E3+18*G2*SCV**3*E1**6-18*E1**6*G2*SCV**2+9*E1**6*G2**2+36*E1**6*SCV**2+18*E1**6*G2*SCV-18*E1**6*SCV*G2**2+9*E1**6*SCV**2*G2**2-18*E1**6*G2')
DISC = F(DISCR)
E2 = (1 + SCV) * E1**2


def map_quantities(mu00, mu11, q01, q10):
    D0 = Matrix(F, [[-mu00 - q01, q01], [q10, -mu11 - q10]])
    D1 = Matrix(F, [[mu00, 0], [0, mu11]])
    A = -D0
    Ai = A.inverse()
    P = Ai * D1
    one = vector(F, [1, 1])
    assert P * one == one, 'P not stochastic'
    # stationary vector of the embedded chain: pi (P - I) = 0, pi.1 = 1
    M = (P.transpose() - Matrix.identity(F, 2))
    M[1, 0] = F(1); M[1, 1] = F(1)
    pi = M.solve_right(vector(F, [0, 1]))
    m = [factorial(k) * (pi * (Ai**k) * one) for k in (1, 2, 3)]
    g2 = P.trace() - 1
    rho1 = ((pi * Ai * P * Ai * one) - m[0]**2) / (m[1] - m[0]**2)
    return m, g2, rho1


def split(expr, tag):
    """expr in F must vanish mod RAD^2-DISC; print A and B of A + RAD*B."""
    num = F(expr).numerator()
    co = num.polynomial(R.gen(0)).list() if num.degree(R.gen(0)) > 0 else [num]
    A = R(0); B = R(0)
    for i, c in enumerate(co):
        c = R(c)
        if i % 2 == 0:
            A += c * DISCR**(i // 2)
        else:
            B += c * DISCR**(i // 2)
    ok = (A == 0 and B == 0)
    alt = ''
    if not ok:
        alt = ' [A^2-B^2*DISC == 0: %s]' % (A**2 - B**2 * DISCR == 0)
    print('  %-34s A==0: %-5s B==0: %-5s -> %s%s'
          % (tag, A == 0, B == 0, 'IDENTITY' if ok else 'NOT ZERO', alt))
    return ok


print('== general branch: rates as rational functions of (E1,SCV,E3,G2,RAD) ==')
mu00 = G2*(-4*E3*G2+4*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3*G2-18*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2-18*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV**2-12*E1**3*G2**2-12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV+12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV+12*E1**3*G2*SCV**2-9*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV+3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)+12*E1**3*G2**2*SCV+9*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2+12*E1**3*G2+12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2-3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**3)/(12*E1**3*G2**3*SCV+3*E1**3*SCV**3*G2-12*E1**3*G2**3+18*E1**3*G2**2*SCV**2-3*E1**3*G2+27*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV**2-9*E1**3*G2*SCV**2+18*E1**3*G2**2-12*E1**3*G2**2*SCV+9*E1**3*G2*SCV-12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**3*SCV-9*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**3*G2-24*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV**2-(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3*SCV**2+4*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3*G2**2+12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**3-(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3+2*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3*SCV+9*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2+24*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV-27*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV+6*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV-12*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2-24*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2+6*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**3-4*E3*G2**2)/E1
mu11 = (-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)
q01 = -3*E1**2*(-6*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV+12*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV-6*G2*SCV*E1**2-3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2+(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3+3*E1**2*G2+6*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2-9*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2*G2+3*E1**2*G2*SCV**2-E3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV-6*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV+6*E1**2*G2**2*SCV+3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2-G2*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3-3*E1**2*G2**2+3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2*G2**2-3*E1**2*SCV**2*G2**2+G2*SCV*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3)/(-45*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV**2+18*G2**2*E1**5*SCV+18*E1**5*G2**3-27*E1**5*G2**2*SCV**2+6*E1**2*G2**2*E3-27*E1**5*G2**2-18*E1**5*G2**3*SCV-18*E1**5*G2*SCV+18*E1**5*G2*SCV**2+3*E1**2*G2*E3-3*E1**2*G2*E3*SCV+(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/E1/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3**2+3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV*E3-36*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV+36*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2+36*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2+45*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV-12*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV*E3-3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*E3+9*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV**3+36*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*SCV**2-6*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**2/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**2*E3+18*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**3*SCV-18*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2**3-9*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)*E1**5/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2)
q10 = 3*(-3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**3-3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV**2+6*E1**3*SCV**2+3*E1**3*G2*SCV**2+3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV**2+6*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2*SCV-E3*SCV-6*E1**3*SCV+(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3*SCV-6*E1**3*G2*SCV-3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*SCV-(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*E3+3*E1**3*G2-3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)*G2+3*E1**3*(-3*E1**3*G2+3*E1**3*G2*SCV-6*E1**3*SCV+E3+RAD)/(-3*E1**3*SCV**2-6*E1**3*SCV-3*E1**3+2*E3)+E3)*E1**2*(-1+G2)/(E3**2-12*E1**3*SCV*E3+6*E1**3*G2*E3-6*G2*SCV*E1**3*E3+18*G2*SCV**3*E1**6-18*E1**6*G2*SCV**2+9*E1**6*G2**2+36*E1**6*SCV**2+18*E1**6*G2*SCV-18*E1**6*SCV*G2**2+9*E1**6*SCV**2*G2**2-18*E1**6*G2)
print('  mu11 =', mu11)
m, g2, rho1 = map_quantities(mu00, mu11, q01, q10)
print('== exact identities modulo RAD^2 - DISC ==')
res = []
res.append(split(m[0] - E1, 'E1(fit) - E1'))
res.append(split(m[1] - E2, 'E2(fit) - (1+SCV)*E1^2'))
res.append(split(m[2] - E3, 'E3(fit) - E3'))
res.append(split(g2 - G2, 'gamma2(fit) - G2'))
res.append(split(rho1 - G2 * (1 - 1 / SCV) / 2, 'rho1(fit) - G2*(1-1/SCV)/2'))
print('  general branch verified:', all(res))

print('== degenerate branch (G2 -> 0, MAP(1) form) ==')
mu00 = 2*(6*E1**3*SCV-E3)/E1/(6*E1**3*SCV+3*E1**3*SCV**2+3*E1**3-2*E3)
mu11 = 0
q01 = 9*E1**5*(SCV-1)*(SCV**2-2*SCV+1)/(6*E1**3*SCV-E3)/(6*E1**3*SCV+3*E1**3*SCV**2+3*E1**3-2*E3)
q10 = -3*(SCV-1)*E1**2/(6*E1**3*SCV-E3)
m, g2, rho1 = map_quantities(mu00, mu11, q01, q10)
r2 = []
r2.append(split(m[0] - E1, 'E1(fit) - E1'))
r2.append(split(m[1] - E2, 'E2(fit) - (1+SCV)*E1^2'))
r2.append(split(m[2] - E3, 'E3(fit) - E3'))
print('  gamma2(fit) =', g2, ' rho1(fit) =', rho1)
print('  degenerate branch verified:', all(r2))

print('== generic MMPP(2) identities (independent of the fit) ==')
S = PolynomialRing(QQ, ['a', 'b', 'p', 'q'], order='lex')
FS = S.fraction_field()
a, b, p, q = [FS(g) for g in S.gens()]
D0 = Matrix(FS, [[-a - p, p], [q, -b - q]])
D1 = Matrix(FS, [[a, 0], [0, b]])
Ai = (-D0).inverse()
P = Ai * D1
one = vector(FS, [1, 1])
M = (P.transpose() - Matrix.identity(FS, 2))
M[1, 0] = FS(1); M[1, 1] = FS(1)
pi = M.solve_right(vector(FS, [0, 1]))
e1 = pi * Ai * one
e2 = 2 * (pi * (Ai**2) * one)
scv = (e2 - e1**2) / e1**2
gg = P.trace() - 1
rr = []
Pk = Matrix.identity(FS, 2)
for k in range(1, 4):
    Pk = Pk * P
    rr.append(((pi * Ai * Pk * Ai * one) - e1**2) / (e2 - e1**2))
print('  g2 = tr(P)-1 =', gg.factor())
print('  SCV - 1      =', (scv - 1).factor())
print('  rho1 - (g2/2)(1-1/SCV) == 0 :', (rr[0] - gg * (1 - 1 / scv) / 2) == 0)
print('  rho2 - g2*rho1 == 0         :', (rr[1] - gg * rr[0]) == 0)
print('  rho3 - g2^2*rho1 == 0       :', (rr[2] - gg**2 * rr[0]) == 0)
idc = scv + gg * (scv - 1) / (1 - gg)
print('  g2 - (IDC-SCV)/(IDC-1) == 0 :', (gg - (idc - scv) / (idc - 1)) == 0)
