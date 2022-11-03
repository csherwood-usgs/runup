# runup_funcs.py

import numpy as np

func_list = ["A17", "D20", "H86", "NH91", "P14", "P16", "R01", "S06", "S11", "V12"]

def calc_L0(Tp = 10.):
    L0 = (9.81*Tp**2)/(2.*np.pi)
    return L0


def A17(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 0.99*Beta*np.sqrt(H0*L0)
    return R2


def D20(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 1.06*(0.0055*np.sqrt(H0*L0)/Beta+(0.32*np.sqrt(H0*L0*Beta)/2.))
    return R2


def H86(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 0.83*Beta*np.sqrt(H0*L0)+0.2*H0
    return R2


def NH91(H0 = 2., L0 = 156.131, Beta = 0.02):
    Hrms = H0/1.416; 
    if Beta >=0.1:
        LR = 0.6*Beta*np.sqrt(Hrms*L0)
        R2 =1.98*LR
    else:
        LR = 0.06*np.sqrt(Hrms*L0)
        R2 =1.98*LR

    return R2


def P14(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 1.29*H0*(Beta/np.sqrt(H0/L0))**(0.72)
    return R2


def P16(H0 = 2., Tp = 10., Beta = 0.02): 
    C = 0.33
    R2 = C*np.sqrt(Beta)*H0*Tp;
    return R2


def R01(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 0.27*np.sqrt(Beta*H0*L0)
    return R2


def S06(H0 = 2., L0 = 156.131, Beta = 0.02):
    Irribian = Beta/np.sqrt(H0/L0)
    S06_Setup = 0.35*Beta*np.sqrt(H0*L0)
    S06_Sinc = 0.75*Beta*np.sqrt(H0*L0)
    S06_Sig = 0.06*np.sqrt(H0*L0)
    S = np.sqrt(S06_Sinc**2+S06_Sig**2)
    R2 = 1.1*(S06_Setup+0.5*S)
    return R2


def S11(H0 = 2., L0 = 156.131, Beta = 0.02):
    R2 = 2.14*np.tanh(0.4*H0);
    return R2


def V12(H0 = 2., L0 = 156.131, Beta = 0.02):
    Irribian = Beta/np.sqrt(H0/L0);
    R2 = 0.53*Beta*np.sqrt(H0*L0)+0.58*Irribian*np.sqrt(H0**3/L0)
    return R2

