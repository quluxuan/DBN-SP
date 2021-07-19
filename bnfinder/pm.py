from numpy import *
import numpy as np
import math
def p_mkde(x, ma, h):
    n = ma.shape[0]
    d = 1
    sxy = np.cov(ma)
    summ = 0
    detS = sxy
    for ix in range(n):
        p20 = np.transpose(x - ma[ix])
        p21 = sxy**(-1)
        p22 = (x - ma[ix])
        p2 = np.dot(p20, p21)
        p2 = np.dot(p2, p22)
        summ = summ + 1 / float(math.sqrt((2 * math.pi) ** d * detS) )* np.exp(-p2 / float(2 * h ** 2))
    pxy = 1 / (n * h ** d) * summ
    return pxy
def p_mkde1(x, ma0, ma1,h):
    ma=np.array([ma0,ma1])
    n = ma0.shape[0]
    d = 2
    sxy = np.cov(np.transpose(ma0),np.transpose(ma1))
    ssxy = mat(sxy)
    summ = 0
    detS = np.linalg.det(sxy)
    for ix in range(n):
        # p20 = np.transpose(x - ma[:, ix])
        # p21 = np.array(linalg.inv(ssxy))
        # p22 = (x - ma[:, ix])
        p2 = np.dot(np.transpose(x - ma[:, ix]), np.array(linalg.inv(ssxy)))
        p2=np.dot(p2,(x - ma[:, ix]))
        #err0=(2 * math.pi) ** d
        summ = summ + 1 / float(math.sqrt(((2 * math.pi) ** d)* detS)) * np.exp(-p2 / float(2 * h ** 2))
    pxy = 1 / float(n * h ** d) * summ
    return pxy

