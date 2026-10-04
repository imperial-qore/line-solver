import numpy as np

def drift_r(dF, ddF, Q, V, W):
    r"""
    Calculates dV and dW following the equations:

    dV[k,s] = \sum dF[k,s,:,:] * V[:,:] + 1/2 \sum ddF[k,s,:,:,:,:] * W[:,:,:,:]
    dW[k,s,k1,s1] = \sum dF[k,s,:,:] * W[:,:,k1,s1] + \sum dF[k1,s1,:,:] * W[:,:,k,s] + Q[k,s,k1,s1]
    """
    dV = np.tensordot(dF, V, axes=([2, 3], [0, 1])) + 0.5 * np.tensordot(ddF, W, axes=([2, 3, 4, 5], [0, 1, 2, 3]))
    dW = np.tensordot(dF, W, axes=([2, 3], [0, 1]))
    dW += np.transpose(dW, axes=[2, 3, 0, 1])
    dW += Q
    return dV, dW


def drift_r_vector(X, N, S, defineDrift, defineDriftDerivativeQ):
    # retrieve mean field (x) and refinement terms (V,W) from state vector
    x = X[0:N*S].reshape((N, S))
    V = X[N*S:2 * (N*S)].reshape((N, S))
    W = X[2 * (N*S):].reshape((N, S, N, S))

    # calculate drift, drift derivatives and Q
    F = defineDrift(x)
    dF, ddF, Q = defineDriftDerivativeQ(x)

    # initialize derivative vector
    dX = np.zeros(2*(N*S) + (N*S)**2)

    # calculate derivatives of V and W
    dV, dW = drift_r(dF, ddF, Q, V, W)

    dX[0:(N*S)] = F.flatten()
    dX[(N*S):2 * (N*S)] = dV.flatten()
    dX[2 * (N*S):] = dW.flatten()
    return dX

