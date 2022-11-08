def calcR2_Raubenheimer(H,T , slope, barslope):
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