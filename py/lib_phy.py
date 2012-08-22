from    math        import log,log10,exp

def calc_Esat(Tair,scheme='goff'):
    EsatWater = {}

    EsatWater['wmo']    = lambda k:10**(
                            10.79574*(1-273.16/k)
                            -5.02800*log10(k/273.16)
                            +1.50475*(10**-4)*(1-10**(-8.2969*(k/273.16-1)))
                            +0.42873*(10**-3)*(10**(4.76955*(1-273.16/k))-1)
                            +0.78614
                            )*100.

    EsatWater['goff']   = lambda k:10**(
                            -7.90298*(373.16/k-1)
                            +5.02808*log10(373.16/k)
                            -1.3816*10**-7*(10**(11.344*(1-k/373.16))-1)
                            +8.1328*10**-3*(10**-3.49149*(373.16/k-1)-1)
                            +log10(1013.246)
                            )*100.

    EsatWater['murp']   = lambda k:exp(
                            54.842763-6763.22/k
                            -4.21*log(k)
                            +0.000367*k
                            +tanh(0.0415*(k-218.8))
                            *(53.878-1331.22/k-9.44523*log(k)+0.014025*k)
                            )

    EsatWater['bolt']   = lambda k:exp(
                            17.67*(k-273.15)/(k+243.5-273.15)
                            )*6.112*100.

    EsatWater['noah']   = lambda k:exp(
                            2501000./461.5*(1./273.15-1./k)
                            )*6.112*100.

    return EsatWater[scheme](Tair)


def calc_Whgt(Wind,obsHgt,outHgt,roughLen):
    '''
    Arya p199??
    '''
    return Wind*log(outHgt/roughLen)/log(obsHgt/roughLen)


def calc_slopeVP(T):
    '''
    * calc. slope vapor pressure curve              [kPa/K]?
    * ref) ??

    T           : observed air temperature (2m by default)  [K]
    '''
    slopeVP = 4098.*(
                     0.6108*exp(17.27*(T-273.12)/T)/
                     T**2.
                     )

    return slopeVP


def calc_psycho(P):
    '''
    * calc. psychometric constant                   [kPa/K]?
    * ref) ??
    '''
    gamma   = 0.665E-3*(P/10.)              # conv. P [hPa] -> [kPa]

    return gamma


def conv_q2e(C,q,P):
    '''
    * conv. specific humidity to vapor pressure
    * ref) ??

    C       : constants
    q       : specific humidity         [??]
    P       : air pressure              [Pa]?
    '''

    w       = q/(1.-q)
    e       = w*P/(C.e-w)

    return e


def r_a(U,C,obsHgtU=2.,obsHgtQ=2.):
    '''
    * calc. aerodynamic resistance
    * ref) http://www.fao.org/docrep/X0490E/x0490e06.htm

      U         : observed wind speed (2m by default)   [m]

      obsHgtU   : height of wind     measurements       [m]
      obsHgtQ   : height of humidity measurements       [m]
    '''

    d       = C.Veg['zeroDis']
    z_om    = C.Veg['roughLenM']
    z_oh    = C.Veg['roughLenH']

    obsHgtU = 2.
    r_a     = log((obsHgtU-d)/z_om) * log((obsHgtQ-d)/z_oh)             \
              /(U*C.k**2)

    return r_a


def r_s(C):
    r_l     = C.Veg['stomataR']
    lai     = C.Veg['lai']
    actLai  = 0.5*lai

    r_s     = r_l/(actLai)
    return r_s


