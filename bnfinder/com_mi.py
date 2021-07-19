import numpy as np
import math
from pm import p_mkde,p_mkde1
def com_mi(X,Y):
    #mutual_information=[]
    X=np.transpose(X)
    Y=np.transpose(Y)
    d = 2
    nx = len(X)
    hx = ((4/(d+2))**(1/float(d+4)))*(nx**(-1/float(d+4)))
    sum1 = 0
    for ii in range(nx):
        xy = np.array([X[ii],Y[ii]])
        xxyy=np.transpose(xy)
        pxy = p_mkde1(xxyy,X,Y,hx)
        px = p_mkde(X[ii],X,hx)
        py = p_mkde(Y[ii],Y,hx)
        sum1 = sum1+math.log(pxy*float(1/(px*py)))
    Ixy = sum1/nx

    #mutual_information.append(Ixy)
    return Ixy

