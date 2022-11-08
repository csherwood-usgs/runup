def wavelength_L(T, h):
    """
    Calculate depth-dependent wavelength
    Input:
       T - wave period [s]
       h - water depth [m]
    Returns:
       L - wavelength [m]
    """
    g = 9.81
    w = T/(2*np.pi)
    kh = qkhfs(w,h)
    L = g*T*T/(2*np.pi)*(np.tanh(kh))
    return L


def iribarren(B, H, T, location="deepwater"):
    """
    Calculate Iribarren wavenumber
    Battjes, 1974
    Description of breaker types
    https://en.wikipedia.org/wiki/Iribarren_number
    Input:
       B - beach slope as fraction rise/run []
       H - wave height at either deepwater or breakpoint [m]
       T - wave period [s]
       location - Either "deepwater" [default] or "breakpoint"
    Returns:
      I - Iribarren number, dimensionless number []
      descr - description of breaker type [string]
    """
    g = 9.81
    Lo = (g*T**2)/(2*np.pi)
    I = B/np.sqrt(H/Lo)

    # description of breaker type
    if location == "deepwater":
        c = np.array((0.5, 3.3))
    else:
        c = np.array((0.4,2.0))

    if I < c[0]:
        descr =  "spilling"
    elif B > c[1]:
        descr = "surging/collapsing"
    else:
        descr =  "plunging"
    return I, descr


def breaker_type(I,location="deepwater"):
    """
    Description of breaker types
    https://en.wikipedia.org/wiki/Iribarren_number
    Input:
       I = Iribarren number, dimensionless number
       location - Either "deepwater" [default] or "breakpoint"
    """
    if location == "deepwater":
        c = np.array((0.5, 3.3))
    else:
        c = np.array((0.4, 2.0))

    if I < c[0]:
        return "spilling"
    elif B > c[1]:
        return "surging/collapsing"
    else:
        return "plunging"


def reverse_shoal(H,T,h):
    """
    Compute offshore wave height by assuming conservation of energy flux E*class
    Nielsen (2009) eqn. 1.7.5
    Input:
        H - Wave height in intermediate depth h (not breaking) [m]
        T - Wave period [s]
        h - Water depth associated with H [m]
    Returns:
        Ho - Wave height in deepwater [m]
    """
    w = 2*np.pi/T
    kh = qkhfs(w,h)
    Ks = 1./( np.sqrt(np.tanh(kh)*(1.+2.*kh/np.sinh(2.*kh))))
    Ho = H/Ks
    return Ho


def jonswap(w: np.ndarray, Hs: float, Tp: float, gamma: float = 3.7) -> np.ndarray:
    '''
    get Jonswap spectra
    :param w: np.ndarray Angular Frequency
    '''
    omega = w
    wp = 2 * np.pi / Tp
    sigma = np.where(omega < wp, 0.07, 0.09)
    a = np.exp(-0.5 * np.power((omega - wp) / (sigma * omega), 2.0))
    sj = 320 * np.power(Hs, 2) * np.power(omega, -5.0) / np.power(Tp, 4) * \
          np.exp(-1950 * np.power(omega, -4) / np.power(Tp, 4)) * np.power(gamma, a)

    return sj


def refract_coeff(theta0, theta1=0.):
    """
    Refraction coefficient from Snells law, assuming straight and parallel isobaths
    Dean and Dalrymple 1984, Eqn. 4.118
    
    Input:
        theta0 - angle of incoming waves relative to coast perpendicular (looking seaward) (radians)
        theta1 - angle of coast perpendicular: zero if wave angle is relative to coast (radian)
        """
    Kr = np.abs( np.cos(theta0)/np.cos(theta1) )
    return Kr


def ursell( Hs, T, h ):
    """
    Calculate Ursell number
    Reussink et al. Eqn 6.
    """
    w = 2*np.pi/T
    kh = qkhfs(w,h)
    k = kh/h
    Ur =0.75*0.5*Hs*k/(kh)**3.
    return Ur


def qkhfs( w, h ):
    """
    Quick iterative calculation of kh in gravity-wave dispersion relationship
    kh = qkhfs(w, h )
    Input
        w - angular wave frequency = 2*pi/T where T = wave period [1/s]
        h - water depth [m]
    Returns
        kh - wavenumber * depth [ ]
    Orbital velocities from kh are accurate to 3e-12 !
    RL Soulsby (2006) \"Simplified calculation of wave orbital velocities\"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
    g = 9.81
    x = w**2.0 *h/g
    y = np.sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = np.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    kh = y
    return kh
