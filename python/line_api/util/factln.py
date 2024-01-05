from scipy.special import gammaln

def factln(x = None):

    y = gammaln(1+x)
    return y
