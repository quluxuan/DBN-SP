import numpy as np
from numpy.linalg import *
import math


def cmi(v1, v2, vcs):
    if vcs:
        v1 = np.matrix(v1)
        v2 = np.matrix(v2)
        vcs = np.matrix(vcs)
        # v1 = v1[:, 0:200]
        # v2 = v2[:, 0:200]
        # vcs = vcs[:, 0:200]
        # v1 = np.transpose(v1)
        # v2 = np.transpose(v2)
        # vcs = np.transpose(vcs)
        v1s = np.vstack((v1, vcs))
        v1scov=np.cov(v1s)
        c1 = np.linalg.det(v1scov)
        v2s = np.vstack((v2, vcs))
        v2scov = np.cov(v2s)
        c2 = np.linalg.det(v2scov)
        vcst = vcs
        [dem1,dem2]=np.shape(vcst)
        if dem1 == 1:
            c3 = np.var(vcst)
        else:
            c3 = np.linalg.det(np.cov(vcst))
        v4s = np.vstack((v1, v2, vcs))
        c4 = np.linalg.det(np.cov(v4s))
        if (c4==0)or(c3==0):
            cmiv = 10
            return cmiv

        else:
            cmiv = 0.5 * math.log((c1 * c2) / (c3 * c4))

        cmiv = abs(cmiv)
        if cmiv == float("inf"):
            cmiv = 1.0e+01
        return cmiv
    else:
        a = np.cov(v1)
        b = np.cov(v2)
        AB = np.cov(v1, v2)
        ab = np.linalg.det(AB)
        mi = 0.5 * math.log(2) * ((a * b) / ab)
        return mi
