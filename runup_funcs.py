# runup_funcs.py

import numpy as np
import os

func_list = ["H86", "N91", "R01", "S06", "S11", "V12", "P14", "A17", "T19", "D20" ]

def which_computer():
    computername = os.environ['COMPUTERNAME']
    drive = 'C:/'
    if computername == 'IGSAGIEGWSCSH10':
        drive = 'D:/'
    return drive, computername


def calc_L0(Tp = 10.):
    L0 = (9.81*Tp**2)/(2.*np.pi)
    return L0


def A17(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Atkinson, A. L., Power, H. E., Moura, T., Hammond, T., Callaghan, D. P., & Baldock, T. E. (2017).
    Assessment of runup predictions by empirical models on non-truncated beaches on the south-east Australian coast.
    Coastal Engineering, 119, 15–31. https://doi.org/10.1016/j.coastaleng.2016.10.001
"""
    R2 = 0.99*Beta*np.sqrt(H0*L0)
    return R2


def D20(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Didier, D., Caulet, C., Bandet, M., Bernatchez, P., Dumont, D., Augereau, E., et al. (2020).
    Wave runup parameterization for sandy, gravel and platform beaches in a fetch-limited, large estuarine system.
    Continental Shelf Research, 192, 104024. https://doi.org/10.1016/j.csr.2019.104024
    eqn. 17
"""
    # R2 = 1.06*(0.0055*np.sqrt(H0*L0)/Beta+(0.32*np.sqrt(H0*L0*Beta)/2.)) # eqn 16 deprecated
    R2 = 0.117*np.sqrt(H0*L0)
    return R2

def T19(H0 = 2., L0 = 156.313, Beta = 0.02):
    """
    Torres-Freyermuth, A., Pintado-Patiño, J. C., Pedrozo-Acuña, A., Puleo, J. A., & Baldock, T. E. (2019).
    Runup uncertainty on planar beaches.
    Ocean Dynamics, 69(11), 1359–1371. https://doi.org/10.1007/s10236-019-01305-y
    eqn. 14
"""
    R2 = 1.76*Beta*np.sqrt(H0*L0)
    return R2


def H86(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Holman, R. A. (1986). Extreme value statistics for wave run-up on a natural beach.
    Coastal Engineering, 9(6), 527–544. https://doi.org/10.1016/0378-3839(86)90002-5
    Table 1b
"""
    Ir = Beta/(np.sqrt(H0*L0))
    R2 = H0*(0.83*Ir + 0.2)
    return R2


def N91(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Nielsen, P., & Hanslow, D. J. (1991).
    Wave Runup Distributions on Natural Beaches.
    Journal of Coastal Research, 7(4), 1139–1152.
    eqns. 9 and 10
    """
    Hrms = H0/1.416; 
    if Beta >=0.1:
        LR = 0.6*Beta*np.sqrt(Hrms*L0)
        R2 =1.98*LR
    else:
        LR = 0.06*np.sqrt(Hrms*L0)
        R2 =1.98*LR

    return R2


def P14(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Paprotny, D., Andrzejewski, P., Terefenko, P., & Furmańczyk, K. (2014). 
    Application of Empirical Wave Run-Up Formulas to the Polish Baltic Sea Coast. 
    PLOS ONE, 9(8), e105437. https://doi.org/10.1371/journal.pone.0105437
    eqn. 6
"""
    R2 = 1.29*H0*(Beta/np.sqrt(H0/L0))**(0.72)
    return R2


def P16(H0 = 2., Tp = 10., Beta = 0.02):
    """
    Poate, T. G., McCall, R. T., & Masselink, G. (2016).
    A new parameterisation for runup on gravel beaches.
    Coastal Engineering, 117, 176–190. https://doi.org/10.1016/j.coastaleng.2016.08.003
    Eqn 9 - gravel beaches (only?)
    """
    C = 0.33
    R2 = C*np.sqrt(Beta)*H0*Tp;
    return R2


def R01(H0 = 2., L0 = 156.131, Beta = 0.02):
    """
    Ruggiero, P., Komar, P. D., McDougal, W. G., Marra, J. J., & Beach, R. A. (2001). 
    Wave Runup, Extreme Water Levels and the Erosion of Properties Backing Beaches. 
    Journal of Coastal Research, 17(2), 407–419.
    Eq. 5
    """
    R2 = 0.27*np.sqrt(Beta*H0*L0)
    return R2


def S06(H0 = 2., L0 = 156.131, Beta = 0.02):
    Ir = Beta/np.sqrt(H0/L0)
    S06_Setup = 0.35*Beta*np.sqrt(H0*L0)
    S06_Sinc = 0.75*Beta*np.sqrt(H0*L0)
    S06_Sig = 0.06*np.sqrt(H0*L0)
    S = np.sqrt(S06_Sinc**2+S06_Sig**2)
    R2 = 1.1*(S06_Setup+0.5*S)
    return R2


def S11(H0 = 2., L0 = 156.131, Beta = 0.02):
    """ 
    Senechal, N., Coco, G., Bryan, K. R., & Holman, R. A. (2011).
    Wave runup during extreme storm conditions.
    Journal of Geophysical Research: Oceans, 116(C7). https://doi.org/10.1029/2010JC006819
    eqn 17.
    """
    R2 = 2.14*np.tanh(0.4*H0)
    return R2


def V12(H0 = 2., L0 = 156.131, Beta = 0.02):
    """Vousdoukas, M. I., Wziatek, D., & Almeida, L. P. (2012). 
    Coastal vulnerability assessment based on video wave run-up observations 
    at a mesotidal, steep-sloped beach. 
    Ocean Dynamics, 62(1), 123–137.
    https://doi.org/10.1007/s10236-011-0480-x
    figure 5a
    """
    Ir = Beta/np.sqrt(H0/L0);
    R2 = 0.53*Beta*np.sqrt(H0*L0)+0.58*Ir*np.sqrt(H0**3/L0)
    return R2


def bar_Ho(T, hb, gamma=0.7 ):
    Hb = hb*gamma
    Ho = reverse_shoal(Hb,T,hb)
    return Ho

def two_slope_IPA_wrapper( Hs, Tp, beta_f, beta_eff):
    """Wrapper for Lange eq9 so you don't have to see spectra
    """
    Ess, f, df = jonswap_ess(Hs, Tp)
    R2 = two_slope_IPA_free(Ess, f, df, beta_f, beta_eff)
    return( R2 )

def jonswap_ess( Hs, Tp, flo = 0.04, fhi = 0.25, npts = 21, acc = 0.05):
    # Generate a Jonswap spectrum in sea-swell band for Lange et al
    # Input:
    #   Hs - Significant wave height (m)
    #   TP - Peak period (s)
    #   flo - Low-frequency cutoff (Hz) [0.04]
    #   fhi - High-frequency cutoff (Hz) [0.25]
    #   npts- number of points [21]
    #   acc - required accuracy [0.02]
    # Returns
    #   Ess - Energy spectrum (m2/Hz)
    #   f   - Frequency (Hz)
    #   df  - delta f (Hz)
    w = np.geomspace(2.*np.pi*flo, 2*np.pi*fhi, npts)
    Ess = 2.*np.pi*jonswap(w, Hs, Tp)
    f = w/(2.*np.pi)
    # check integration
    dff = np.diff(f)
    df = np.append(dff[0], dff)
    Hsf = 4.*np.sqrt( np.sum(Ess*df))
    err = np.abs(Hsf-Hs)/Hs
    assert err <= acc,'Hs does not match input'
    return Ess, f, df


def two_slope_IPA_free(Ess, f, df, beta_f, beta_eff):
    """
    Equation 9 in Lange et al., 2020
    """
    #term1 = 0.4 * 𝛽eff^0.45 ∫𝑆𝑆 𝐸^0.5 𝑓^−1.45 𝑑𝑓
    aint = Ess**0.5 * f**-1.45 * df
    term1 = 0.4 * beta_eff**0.45 * aint.sum()

    #term2 = 3.72 * 𝛽eff^0.6 ∫𝑆𝑆 𝐸^0.95 𝑓^−0.25 𝑑𝑓
    bint = Ess**0.95 * f**-0.25 * df
    term2= 3.72 * beta_eff**0.6 * bint.sum()

    #term3 = 0.76 * 𝛽f^2 * 𝛽eff^0.75 ∫𝑆𝑆 𝐸^0.4 𝑓^−3.1 𝑑𝑓
    cint = Ess**0.4 * f**-3.1 * df
    term3= 0.76 * beta_f**2 * beta_eff**0.75 * cint.sum()

    #R2%,G = term1 + 2[term2 + term3]^0.5
    R2 = term1 + 2*(term2+term3)**0.5
    return R2


def ang_corr (H, a, an):
    """Correct for angle of approach
    Input: 
        H = wave height
        a = approach angle (degrees)
        an = shore normal angle (degrees)
    Returns:
        Hc = corrected wave height
    """
    adiff = np.deg2rad( np.abs(a-an) )
    Hc = np.sqrt( np.cos (adiff) *H*H )
    return Hc

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
        % Converted to Python by csherwood@usgs.gov

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


def calcR2(H, T, slope, igflag=0):
    """
    %
    % [R2, S, setup,  Sinc,  SIG,  ir] = calcR2(H, T, slope, igflag);
    %
    % Calculated 2% runup (R2), swash (S), setup (setup), incident swash (Sinc)
    % and infragravity swash (SIG) elevations based on parameterizations from runup paper
    % also Iribarren (ir)
    % August 2010 - Included 15% runup (R16) statistic that, for a Guassian distribution,
    % represents mean+sigma. It is calculated as R16 = setup + swash/4.
    % In a wave tank, Palmsten et al (2010) found this statistic represented initiation of dune erosion.
    %
    %
    % H = significant wave height, reverse shoaled to deep water
    % T = deep-water peak wave period
    % slope = radians
    % igflag = 0 (default)use full equation for all data
    %        = 1  use dissipative-specific calculations when dissipative conditions exist (Iribarren < 0.3)
    %        = 2  use dissipative-specific (IG energy) calculation for all data
    %
    % based on:
    %  Stockdon, H. F., R. A. Holman, P. A. Howd, and J. Sallenger A. H. (2006),
    %    Empirical parameterization of setup, swash, and runup,
    %    Coastal Engineering, 53, 573-588.
    % author: hstockdon@usgs.gov
    # Converted to Python by csherwood@usgs.gov
    """
    g = 9.81

    # make slopes positive!
    slope = np.abs(slope)

    # compute wavelength and Iribarren
    L = (g*T**2) / (2.*np.pi)
    sqHL = np.sqrt(H*L)
    ir = slope/sqHL

    if igflag == 2:                     # use dissipative equations (IG) for ALL data
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    elif igflag == 1 and ir < 0.3:      # if dissipative site use diss equations
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    else:                               # if int/ref site, use full equations
        setup = 0.35*slope*sqHL
        Sinc = 0.75*slope*sqHL
        SIG = 0.06*sqHL
        S = np.sqrt(Sinc**2 + SIG**2)
        R2 = 1.1*(setup + S/2.)
        R16 = 1.1*(setup + S/4.)

    return R2, S, setup, Sinc, SIG, ir, R16


def stat_summary(x, iprint=False):
    n = len(x)
    nnan = np.sum(np.isnan(x))
    nvalid = n-nnan
    # intitialize with NaNs

    if n > nnan:
        meanx = np.nanmean(x)
        stdx = np.nanstd(x)
        minx = np.nanmin(x)
        d5 = np.nanpercentile(x, 5.)
        d25 = np.nanpercentile(x, 25.)
        d50 = np.nanpercentile(x, 50.)
        d75 = np.nanpercentile(x, 75.)
        d95 = np.nanpercentile(x, 95.)
        maxx = np.nanmax(x)
    else:
        meanx = np.nan
        stdx = np.nan
        minx = np.nan
        d5 = np.nan
        d25 = np.nan
        d50 = np.nan
        d75 = np.nan
        d95 = np.nan
        maxx = np.nan

    # return it in a dict
    s = {'n':n, 'nnan':nnan, 'nvalid':nvalid, 'mean':meanx, 'std':stdx, 'min':minx, 'max':maxx,
         'd5':d5, 'd25':d25, 'd50':d50, 'd75':d75, 'd95':d95}
    # if iprint:
    #     for key, value in s.items():
    #         print('{:6s} = {:.3f}'.format(key, value)),
    if iprint:
        print("  n, nnan, nvalid: ",s['n'],s['nnan'],s['nvalid'])
        print("  mean, std, min, max   : {:.3f} {:.3f} {:.3f} {:.3f}"
            .format(s['mean'], s['std'], s['min'], s['max']))
        print("  d5, d25, d50, d75, d95: {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}"
            .format(s['d5'], s['d25'], s['d50'], s['d75'], s['d95']))

    return s
