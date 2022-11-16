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


def calcR2_Raubenheimer(H, T , slope, barslope):
        """
        %
        % [R2R,S,setup, Sinc, SIG, ir, R16R] = calcR2_Raubenheimer(H,T,slope,barslope);
        %
        % Calculated 2% runup (R2R), swash (S), setup (setup), incident swash (Sinc)
        % and infragravity swash (SIG) elevations based on parameterizations from runup paper
        % also Iribarren (ir)
        % August 2010 - Included 15% runup (R16R) statistic that, for a Gaussian distribution, 
        % represents mean+sigma. It is calculated as R16 = setup + swash/4.  
        % In a wave tank, Palmsten et al (2010) found this statistic represented initiation of dune erosion. 
        %
        %
        % R2R=1.3 * (2.1*barslope(H0*L0)^1/2+ [(0.005 * slope^2 * L0^2 + 0.004*H0*L0)]^1/2 ?2)
        %
        %
        % H = significant wave height, reverse shoaled to deep water
        % T = deep-water peak wave period
        % slope = radians
        % barslope = radians
        %
        % based on:
        % Britt Raubenheimer, personal communication
        % From the Matlab version provided by Alfredo
        """
        g = 9.81

        # make slopes positive
        slope = np.abs(slope)
        barslope = np.abs(barslope)

        # compute wavelength and Iribarren
        L = (g*T**2) / (2.*np.pi)
        ir = slope/np.sqrt(H/L)

        setup = 2.1*barslope*np.sqrt(H*L)
        Sinc = np.sqrt(0.005)*slope*L
        SIG = np.sqrt(0.004*H*L)
        S = np.sqrt(Sinc**2 + SIG**2)
        R2R = 1.3*(setup + S/2.)
        R16R = 1.3*(setup + S/4.)
        return R2R, S, setup, Sinc, SIG, ir, R16R