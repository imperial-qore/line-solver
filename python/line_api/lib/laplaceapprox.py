def laplaceapprox(h = None,x0 = None):
    # I = laplaceapprox(f,x0)
    # approximates I=int f(x)dx by Laplace approximation at x0
    # example:  I = laplaceapprox(@(x) prod(x),[0.5,0.5])
    d = len(x0)

    tol = 1e-05
    H = num_hess(lambda x = None: log(h(x)),x0,tol)
    detnH = det(- H)
    if detnH < 0:
        tol = 0.0001
        H = num_hess(lambda x = None: log(h(x)),x0,tol)
        detnH = det(- H)

    if detnH < 0:
        tol = 0.001
        H = num_hess(lambda x = None: log(h(x)),x0,tol)
        detnH = det(- H)

    if detnH < 0:
        warnings.warn('laplaceapprox.m: det(-H)<0')

    I = h(x0) * sqrt((2 * pi) ** d / detnH)
    logI = log(h(x0)) + (d / 2) * log(2 * pi) - log(detnH)
    return I,H,logI
