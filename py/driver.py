#! /usr/bin/python
#--------------------------------------
# PROG : driver.py	 <by hjkim@IIS>
# VER  : <2012-07-21 12:30:51.424602>
#
# DESC :
# USAGE: $ ./driver.py 
#------------------------cf@none


import  os,sys,datetime
from    math    import log10,log

from    const   import Const
from    lib_phy import calc_Esat,calc_Whgt,calc_slopeVP,calc_psycho,    \
                       conv_q2e,r_a,r_s

#from pylab import *



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


def calc_H(C,dVarIn,dVarState):
    '''
    * calc. sensible heat flux
    * ref) 

    H               : Sensible Heat [W/m^2]

    C               : Constants
    dVarIn          : Input Variables
    '''

    Cp      = C.Cp
    Cd      = C.Cd

    PSurf   = dVarIn['PSurf']*100.      # unit conv. [hPa] -> [Pa]
    Tair    = dVarIn['Tair']
    U10     = dVarIn['Wind']

    rho     = C.rho(PSurf,Tair)

    Ts      = dVarState['Ts']

    H       = Cp*rho*Cd*U10*(Ts-Tair)

    return H


def calc_ET0(C,slopeVP,Rnet,G,gamma,T2,U2,Q2,P):
    '''
    * calc. reference evapotranspiration using FAO P-M eq.
    * ref) FAO ch2

    slopeVP : slope vapor pressure curve        [kPa/K]
    Rnet    : net radiation                     [W/m**2]
    G       : soil heat flux                    [W/m**2]
    gamma   : psychometric constant             [kPa/K]
    T2      : 2m air temperature                [K]
    U2      : 2m wind speed                     [m/s]
    Q2      : 2m specific humidity              [??]
    P       : surface pressure                  [hPa]
    '''

    Esat    = calc_Esat(T2,scheme='goff')
    E       = conv_q2e(C,Q2,P*100.)         # conv. unit [hPa] -> [Pa]

    # VPD     : saturation vapor pressure deficit [kPa]
    VPD     = (Esat-E)/1000.                # conv. unit [hPa] -> [kPa]

    Rnet    = Rnet*86400/10**6              # conv. unit [W/m^2] -> [MJ/m^2/d]
    G       = G*86400/10**6                 # conv. unit [W/m^2] -> [MJ/m^2/d]

    ET0     = (0.408*slopeVP*(Rnet-G)+gamma*900./(T2)*U2*VPD)/  \
              (slopeVP+gamma*(1.+0.34*U2))

#    print ('%g '*8)%(slopeVP,Rnet,G,gamma,T2,U2,Q2,P);print ET0;print

    return ET0


#def update_Ts(C.Veg['albedo'],RSDN,RLDN,C.sig,ET0,H):
def update_Ts(C,dVarIn,ET0,H):

    albedo  = C.Veg['albedo']
    
    RSDN    = dVarIn['SWdown']
    RLDN    = dVarIn['LWdown']

    sigma   = C.sig

    print RSDN,RLDN
    print ET0, H
    print (1.-albedo)*RSDN
    print sigma
    print ((1.-albedo)*RSDN+RLDN-ET0-H)
    Ts      = (((1.-albedo)*RSDN+RLDN-ET0-H)/sigma)**(1./4.)

    return Ts


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

    # Open Input Files --------------------------------------------------------
    dInFile     = dict(
                    (var,file(os.path.join(
                                inputDir,
                                '%s.%s.%i@%ix%i.asc'%(prjName,var,sDTime.year,yIdx,xIdx)
                                           ))) 
                                for var in varNAME)
    # --------------------------------------------------------------------------
    # Declare State Variables
    dVarState   = {
                   'Ts':273.15, # Skin Temp.
                   'Td':274.15, # Soil Temp.
                    }


    tmpOUT  = {'ET0':[]}
    # Time integration loop ---------------------------------------------------
    for nLoop   in xrange(nTLoop):
        dVarIn  = dict((var,float(dInFile[var].readline())) 
                                for var in varNAME)
        
#        Epot    = calc_Epot(C,dVarIn)
#        print '%5.2f'%(Epot*86400), 

        # Read State Variable -------------------------------------------------
        Ts      = dVarState['Ts']               # Skin Temp.            [K]
        # ---------------------------------------------------------------------

        # Read Forcing Variable -----------------------------------------------
        U10     = dVarIn['Wind']                # 10m wind speed        [m/s]
        T2      = dVarIn['Tair']                # 2m  air temp.         [K]
        Q2      = dVarIn['Qair']                # 2m specific humidity  [??]
        P       = dVarIn['PSurf']               # surface pressure      [hPa]

        RSDN    = dVarIn['SWdown']              # downward solar rad.   [W/m**2]
        RLDN    = dVarIn['LWdown']              # downward solar rad.   [W/m**2]
        # ---------------------------------------------------------------------

        RSUP    = RSDN*C.Veg['albedo']          # upward   solar rad.   [W/m**2]
        RLUP    = C.sig*Ts**4                   # upward   solar rad.   [W/m**2]

        Rnet    = RSDN-RSUP+RLDN-RLUP

        U2      = calc_Whgt(U10,10.,2.,C.Veg['roughLenM'])    # calc. 2m wind speed

        slpVP   = calc_slopeVP(T2)
        gamma   = calc_psycho(P)

        R_a     = r_a(U2,C,obsHgtU=2.)
        R_s     = r_s(C)

#        print dVarIn['Wind'],U10,r_a(U2,C,obsHgtU=2.),207.66407000788683/U2,69.44444444
#        print T2,U10,R_a,R_s,RSDN,RSUP,RLDN,RLUP,Rnet,gamma

        ET0     = calc_ET0(C,slpVP,Rnet,0.,gamma,T2,U2,Q2,P)
        H       = calc_H(C,dVarIn,dVarState)

        dVarState['Ts'] = update_Ts(C,dVarIn,ET0,H)

        tmpOUT['ET0'].append(ET0)

    return  tmpOUT


if __name__=='__main__':
    tmpOUT  = main(*sys.argv)

    from pylab import *
    aTmp    = array(tmpOUT['ET0']).reshape(-1,4).mean(1)
    plot(aTmp);show()
        

    
