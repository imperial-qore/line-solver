"""
Fold auxiliary-class metrics back into original classes (native fork-join).

Port of matlab/src/api/fj/sn_fj_foldback.m.
"""

import numpy as np


def sn_fj_foldback(QN, UN, RN, TN, CN, XN, fjclassmap, Korig):
    """
    Fold the auxiliary-class columns of the average metrics computed on an
    FJ tag-augmented struct (ModelAdapter.fjtag) back into the original
    classes: queue lengths, utilizations and throughputs of the sibling
    classes are exact aggregates of the original class they were forked
    from; response times are recomputed by Little's law after folding.

    QN,UN,RN,TN are (nstations x Kaug); CN,XN are (1 x Kaug); the outputs
    retain only the first Korig columns.
    """
    QN = np.array(QN, dtype=float, copy=True)
    UN = np.array(UN, dtype=float, copy=True)
    RN = np.array(RN, dtype=float, copy=True)
    TN = np.array(TN, dtype=float, copy=True)
    CN = np.atleast_2d(np.array(CN, dtype=float, copy=True))
    XN = np.atleast_2d(np.array(XN, dtype=float, copy=True))

    fjclassmap = np.asarray(fjclassmap).ravel()
    Kaug = len(fjclassmap)
    for a in range(Kaug):
        r = int(fjclassmap[a])
        # fjclassmap stores the 0-based original class of each auxiliary
        # class, or -1 for original classes (python 0-based convention;
        # MATLAB uses 0 as the "original" sentinel over 1-based indices).
        if r >= 0:
            QN[:, r] = QN[:, r] + QN[:, a]
            UN[:, r] = UN[:, r] + UN[:, a]
            TN[:, r] = TN[:, r] + TN[:, a]

    QN = QN[:, :Korig]
    UN = UN[:, :Korig]
    TN = TN[:, :Korig]
    RN = np.zeros_like(QN)
    for r in range(Korig):
        for ist in range(QN.shape[0]):
            if TN[ist, r] > 0:
                RN[ist, r] = QN[ist, r] / TN[ist, r]

    # system metrics are measured on the original class columns only
    CN = CN[:, :Korig]
    XN = XN[:, :Korig]

    return QN, UN, RN, TN, CN, XN
