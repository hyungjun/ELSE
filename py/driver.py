#! /usr/local/bin/python
#--------------------------------------
# PROG : driver.py	 <by hjkim@IIS>
# VER  : <2012-07-21 12:30:51.424602>
#
# DESC :
# USAGE: $ ./driver.py 
#------------------------cf@none


import  os,sys,datetime
from    math    import log10,log

from    lib_phy import calc_Esat, calc_Whgt

#from pylab import *

class Const(object):    # constants

    Cd      = 0.003     # bulk transfer coefficient
    e       = 0.62198   # ratio of molecular weights of water and dry air
    Cp      = 10005     # Specific Heat Capacity of Air [J/Kg/K]
    k       = 0.41      # von Karman constant [-]

        
    dictVeg = {'crop':{'height'   :0.12,          # veg. height [m]
                       'zeroDis'  :2./3.*0.12,    # zero plane displacement              [m] 2/3*H
                       'roughLenM':0.123*0.12,    # roughness length for momentum   flux [m] 0.123*H
                       'roughLenH':0.123*0.12*0.1,# roughness length for heat&vapor flux [m] 0.123*H*0.1
                      }
              }

    def __init__(self,vegType='crop'):
        self.Veg        = self.dictVeg[vegType]

        print self.Veg

    # calc. air density 
    def rho(self,PSurf=101325.,Tair=273.12+15.):    
        '''
        rho     = 1.225     # air density       [kg/m^3] (sea level 15degC by ISA)

        PSurf               # air pressure      [pa]
        Tair                # air temperature   [K]

        '''
        Rg      = 287.04    # specific gas constant [J/kg/K]
        return PSurf/Tair/Rg


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


def calc_Epot(C,dVarIn,Esat_scheme='goff'):
    '''
    Epot            : potential evaporation [kg/m^2/s]

    C               : Constants
    dVarIn          : Input Variables
    '''

    e       = C.e
    Cd      = C.Cd

    PSurf   = dVarIn['PSurf']*100.      # unit conv. [hPa] -> [Pa]
    Tair    = dVarIn['Tair']
    U       = dVarIn['Wind']
    Qair    = dVarIn['Qair']

    rho     = C.rho(PSurf,Tair)
    Esat    = calc_Esat(Tair,scheme=Esat_scheme)
    Qsat    = e*Esat/(PSurf-Esat)

    Epot    = rho*Cd*U*(Qsat-Qair)

    return Epot


def calc_H(C,dVarIn):
    '''
    H               : Sensible Heat [W/m^2]

    C               : Constants
    dVarIn          : Input Variables
    '''

    Cp      = C.Cp
    Cd      = C.Cd

    PSurf   = dVarIn['PSurf']*100.      # unit conv. [hPa] -> [Pa]
    Tair    = dVarIn['Tair']
    U       = dVarIn['Wind']
    Qair    = dVarIn['Qair']

    rho     = C.rho(PSurf,Tair)


def main(*args):

    inputDir    = './data.1d'

    prjName     = 'Prcp_GPCC'
    xIdx,yIdx   = 100,140

    vegType     = 'crop'

    C           = Const(vegType)

    sDTime      = datetime.datetime(2000,1,1,0,0)
    eDTime      = datetime.datetime(2001,1,1,0,0)
    dT          = datetime.timedelta(seconds=3600*6)

#    totSec      = (eDTime-sDTime).days*86400+(eDTime-sDTime).seconds
#    dTsec       = dT.days*86400+dT.seconds
    totSec      = (eDTime-sDTime).total_seconds()   # in float (not int)
    dTsec       = dT.total_seconds()

    nTLoop      = int(totSec/dTsec)

    varNAME     = [
                   'CCOV',
                   'LWdown',
                   'PSurf',
                   'Prcpf',
                   'Qair',      # 2m Specific Humidity  [kg/kg]
                   'Rainf',
                   'SWdown',
                   'Snowf',     
                   'Tair',      # 2m Air Temperature    [K]
                   'Wind'       # Wind Speed            [m/s]
                   ]

    dInFile     = dict(
                    (var,file(os.path.join(
                                inputDir,
                                '%s.%s.%i@%ix%i.asc'%(prjName,var,sDTime.year,yIdx,xIdx)
                                           ))) 
                                for var in varNAME)


    for nLoop   in xrange(nTLoop):
        dVarIn  = dict((var,float(dInFile[var].readline())) 
                                for var in varNAME)
        
#        Epot    = calc_Epot(C,dVarIn)
#        print '%5.2f'%(Epot*86400), 

        U10     = dVarIn['Wind']
        U2      = calc_Whgt(U10,10.,2.,C.Veg['roughLenM'])    # calc. 2m wind speed

        print dVarIn['Wind'],U10,r_a(U2,C,obsHgtU=2.),207.66407000788683/U2
#        print dVarIn['Wind'],U,r_a(U2,C,obsHgtU=2.),208./U2

    return


if __name__=='__main__':
    main(*sys.argv)
    print 'test'

    
