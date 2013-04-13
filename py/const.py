from    numpy   import sqrt


class Const(object):        # constants

    Cd      = 0.003         # bulk transfer coefficient         ref) Robock et al., 1995?
    e       = 0.62198       # ratio of molecular weights of water and dry air
    Cp      = 1005          # Specific Heat Capacity of Air [J/Kg/K]
    k       = 0.41          # von Karman constant [-]
    sig     = 5.670373E-8   # Stefan-Boltzmann constant [W/m**2/K**4]

    Cs      = 13000.        # heat capacity of the surface layer [J/m**2/K]
    Cd_     = sqrt(365.)*Cs # heat capacity of the deep layer [J/m**2/K]
    w       = 1./(24.*60.*60.)  # time constant [1/s]

        
    dictVeg = {'crop':{'height'   :0.12,          # veg. height [m]
                       'zeroDis'  :2./3.*0.12,    # zero plane displacement              [m] 2/3*H
                       'roughLenM':0.123*0.12,    # roughness length for momentum   flux [m] 0.123*H
                       'roughLenH':0.123*0.12*0.1,# roughness length for heat&vapor flux [m] 0.123*H*0.1
                       'stomataR' :100.,          # well-watered stomatal resistance     [s/m]
                       'lai'      :24*0.12,       # leaf area index                      [-] 24*H
                       'albedo'   :0.23           # albedo
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

