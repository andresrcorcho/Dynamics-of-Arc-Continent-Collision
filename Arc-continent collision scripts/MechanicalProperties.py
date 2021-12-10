import underworld as uw
import math
from underworld import function as fn
import numpy as np
import scipy
import os
import matplotlib.pyplot as pypl
import scipy.special as special
import UWGeodynamics as GEO

# (2.6+5)/2.
# 3.8


# +
#3546.2

# +
#3400+(3400*0.043)
#Mean between the limits 1.6 and 7, (including errors)

# +
#3400+(3400*0.038)
#mean between the limits 2.6 and 5 %
# -

# 3400+(3400*0.026)
# 3488.4 - dwnsity of lower mantle with lowest contrast of 2.6 %


u = GEO.UnitRegistry

#Model Parameters
Tsurf = 273.15 * u.degK
Tint  = 1573.0 * u.degK
kappa = 1e-6   * u.meter**2 / u.second 
alpha = 3.0e-5 / u.kelvin

grav  = 9.81   * u.meter / u.second**2
R     = 8.314  * u.joule / u.degK / u.mol

# wet diffusion creep - Karato & Wu, 1993
A = 5.3  * 10**15 / u.pascal / u.second
E = 240. * 10**3 * u.joule / u.mol
V = 5.0  * 1e-6 * u.meter**3 / u.mol

max_viscosity = 10**5

def density_defArc(depth,crustDensity,mantleDensity,crustThickness,rho,arcDensity,arcThickness):
    for index, y in enumerate(depth):
        if y < crustThickness:
            rho[index] = crustDensity
        elif y < arcThickness:
            rho[index] = arcDensity
        else:
            rho[index] = mantleDensity
    return rho

def density_def(depth,crustDensity,mantleDensity,crustThickness,rho,arc,arcDensity,arcThickness):
    if arc==True:
        rho=density_defArc(depth,crustDensity,mantleDensity,crustThickness,rho,arcDensity,arcThickness)
    else:
        for index, y in enumerate(depth):
            if y < crustThickness:
                rho[index] = crustDensity
            else:
                rho[index] = mantleDensity
    return rho

#Arrhenius Viscosity
def arrhenius(T,P):
    """ Arrhenius viscosity assuming absolute temperature T and pressure P """
    #print( E, V, P, R, T)
    # DW should be scaled to 'density' of crust or mantle, not just mantle?
    return np.exp((E + V*P)/(R*T))

# for oceanic lithosphere thickness
def half_space_cooling(depth,age):
    return Tsurf +(Tint-Tsurf) * special.erf(depth/(2*np.sqrt(age*kappa)))

# temperature of continental lithosphere
def linear_geotherm(depth, plateThickness):
    linear = depth.to('kilometer')*(Tint-Tsurf)/(plateThickness*u.kilometer)+Tsurf
    return np.minimum(linear.magnitude,Tint)*linear.units

def yeildStrength(depth,crustThickness,plateThickness,cohesionC, cohesionLit,friction,friction2,pressure,cohesionFactor):
    yieldStrength = np.zeros(len(depth))*u.megapascal
    yieldStrength2 = np.zeros(len(depth))*u.megapascal
    for index, y in enumerate(depth):
        if y < crustThickness: #In Kilometers
            yieldStrength[index] =  cohesionC + friction  * pressure[index]
            yieldStrength2[index] =  (cohesionC/cohesionFactor) + (friction2)  * pressure[index] #Post-yield strength
        elif y < plateThickness:  #Varying this makes the overriding plate weak
            yieldStrength[index] =  cohesionLit + friction  * pressure[index]
            yieldStrength2[index] =  (cohesionLit/cohesionFactor) + (friction2) * pressure[index] #Post-yield strength
    return yieldStrength,yieldStrength2

def PlateProperties(Nlayers,crustThickness,crustDensity,mantleDensity,plateThickness,oceanic,age, cohesionC,cohesionLit,friction,friction2,arc,arcDensity,arcThickness,depthToMantle,cohesionFactor):
    #Model Parameters
    age=age * 10**6 * u.year #Megayears
    res   = 251
    depth = np.linspace(0, 250., res) * 10**3 * u.meter
    mDep=250.
    prop=plateThickness/Nlayers
    layers= np.array([])
    depAc=0.
    for i in range(0,Nlayers+1):
        layers=np.append(layers,[depAc])
        depAc=depAc+prop
    layers=np.append(layers,[mDep])* 10**3 * u.meter 
    #Densities in depth and crust thc
    crustThickness=crustThickness* u.kilometer
    plateThickness=plateThickness* u.kilometer
    crustDensity=crustDensity  * u.kilogram / u.meter**3
    mantleDensity=mantleDensity  * u.kilogram / u.meter**3
    arcDensity=arcDensity * u.kilogram / u.meter**3
    arcThickness=arcThickness * u.kilometer
    #Reference Density
    rho = np.ones(res)* u.kilogram / u.meter**3
    #Density profile
    rho=density_def(depth,crustDensity,mantleDensity,crustThickness,rho,arc,arcDensity,arcThickness)
    #Density
    if oceanic==True:
        temperature = half_space_cooling(depth.to_base_units(), age.to_base_units())
    else:
        temperature = linear_geotherm(depth,depthToMantle)
        
    deltaT  = temperature - Tsurf
    density = rho - rho * deltaT * alpha
    # lithostatic pressure = \rho * g * h
    pressure = density * grav * depth
    #Mechanical strength
    cohesionC = cohesionC * u.megapascal
    cohesionLit=cohesionLit* u.megapascal
    
    yieldStrength  = yeildStrength(depth,crustThickness,plateThickness,cohesionC, cohesionLit,friction,friction2,pressure,cohesionFactor)[0]
    yieldStrength2 = yeildStrength(depth,crustThickness,plateThickness,cohesionC, cohesionLit,friction,friction2,pressure,cohesionFactor)[1]
    
    # define some reference parameters
    ref_depth     = 100. * 10**3 *u.meter   
    ref_pressure  = mantleDensity * grav * ref_depth
    ref_temp      = Tint
    ref_viscosity = arrhenius(ref_temp, ref_pressure)
    ref_density   = mantleDensity * (Tint-Tsurf) * alpha
    # calculate viscosity profile
    viscosity = np.clip(arrhenius(temperature, pressure)/ref_viscosity, a_min = 0, a_max = max_viscosity) 
    # calculate averages of layers
    avg_depth      = []
    avg_temp       = []
    avg_density    = []
    avg_pressure   = []
    avg_viscosity  = []
    avg_strength   = []
    low_strength   = []
    avg_strength2  = []
    low_strength2  = []

    prevIndex = 0
    for i, interface in enumerate(layers):
        if( interface > 250.*u.kilometer):
            break
        index = np.where(depth.magnitude==interface.magnitude)[0][0]
        if depth[index].magnitude > 0.:
            avg_depth.append(    np.average(depth.to('km').magnitude[prevIndex:index]))
            avg_temp.append(     np.average(temperature[prevIndex:index]))
            avg_density.append(( np.average(density[ prevIndex: index])-(density[-1]).magnitude)/ref_density.magnitude)
            avg_pressure.append( np.average(pressure.to('MPa').magnitude[prevIndex:index]))
            avg_viscosity.append(np.average(viscosity.magnitude[prevIndex:index]))
            avg_strength.append( np.average(yieldStrength.to('MPa').magnitude[prevIndex:index]))
            low_strength.append( np.min(yieldStrength.to('MPa').magnitude[prevIndex:index]))
            avg_strength2.append( np.average(yieldStrength2.to('MPa').magnitude[prevIndex:index]))
            low_strength2.append( np.min(yieldStrength2.to('MPa').magnitude[prevIndex:index]))

        prevIndex = index
        
        densities=[]
        viscosities=[]
        for i in avg_density:
            densities.append((i*ref_density.magnitude)+(mantleDensity.magnitude))
        for i,vis in enumerate (avg_viscosity):
            viscosities.append(vis*1e20)
        #return avg_density, avg_viscosity, avg_strength,avg_strength2, layers
        
    return avg_density, avg_viscosity, densities,viscosities,avg_strength,avg_strength2,low_strength,low_strength2,[depth, temperature,viscosity,layers,density,age,pressure,yieldStrength,yieldStrength2]
