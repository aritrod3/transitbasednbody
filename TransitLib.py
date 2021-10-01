import numpy as np
import matplotlib.pyplot as plt
from Particle import Particle
import pdb
import astropy.units as u
plt.ion()
import pyvo as vo
from astropy.table import QTable
import astropy.units as u


def inputPlanet():
    NArx_Service = vo.dal.TAPService('https://exoplanetarchive.ipac.caltech.edu/TAP')
    col2pull =  "pl_name,hostname,pl_orbsmax,pl_orbeccen,pl_orbper,pl_bmasse," + \
                    "st_mass," + \
                    "discoverymethod"

    tab2pull = 'pscomppars'
    query = 'select '+col2pull+' from '+tab2pull
    NArx_res = NArx_Service.search(query)
    qtab = (QTable(NArx_res.to_table()))[(QTable(NArx_res.to_table()))['discoverymethod'] == 'Transit']
    
    OrbitTable = qtab[qtab['pl_orbeccen'] == 0.00]
    
    print(OrbitTable)
    exo = input("Give a transiting exoplanet:")
    
    OrbitTable = OrbitTable[OrbitTable['pl_name'] == exo]
    
    
    col2pull = 'plntname,centralwavelng,bandwidth,plntransdep,plnradj'
    tab2pull = 'transitspec'
    TransitData = NArx_Service.search(query)
    qtab = (QTable(NArx_res.to_table()))[(QTable(NArx_res.to_table()))['discoverymethod'] == 'Transit']
    TransitTable = QTable(TransitData.to_table())
    
    
    OrbitTable = OrbitTable[OrbitTable['pl_name'] == exo]
    print(OrbitTable)
    print("---")
    print(OrbitTable['pl_orbeccen'])
    
    
    if(exo not in OrbitTable['pl_name']):
        raise KeyError("Exoplanet not found")
    else:
        StellarMass = ((OrbitTable['st_mass']).si.value)[0]
        print("Stellar Mass = " + str(StellarMass) + " kg")
        
        PlanetMass = (((OrbitTable['pl_bmasse']).value)[0])*(5.972e24)
        print("Planetary Mass = " + str(PlanetMass) + " kg")
        PlanetarySemiMajorAxis = (((OrbitTable['pl_orbsmax']).to(u.m).value)[0])
        print("Planetary Semi Major Axis  = " + str(PlanetarySemiMajorAxis) + " m")
        PlanetOrbPer = (((OrbitTable['pl_orbper']).value)[0]*(86400))
        print("Planetary Orbital Period = " + str(PlanetOrbPer) + " s")
        PlanetVelo = 1.4*2*np.pi*PlanetarySemiMajorAxis/PlanetOrbPer
        print("Planetary Orbital Velocity = " + str(PlanetVelo) + " m/s")
    
    return (StellarMass, PlanetMass, PlanetarySemiMajorAxis, PlanetVelo, PlanetOrbPer)
        
#Normalizes the flux between 0 and 1, where 1 is the brightness of the original star
def normalize(arr, star_matrix):
    for i in range(len(arr)):
        arr[i] /= (np.sum(star_matrix))
    return arr


#This updates the transmittance such that the runTransit method can calculate the fluxes in different atmospheric wavelengths.
def updateTransmittance(planet_matrix, new_transmittance):
    u = np.unique(planet_matrix)
    t = u[-2]
    
    planet_matrix[planet_matrix==t] = new_transmittance
    
    return planet_matrix
                
                


                
#This runs the transit
def runTransit(hShift, hIncr, vShift, vIncr, steps, star_matrix, planet_matrix, transmissions):

    for i in range(steps):
        if hIncr != 0 : planet_matrix = np.roll(planet_matrix, hIncr, 1)
        if vIncr != 0: planet_matrix = np.roll(planet_matrix, -vIncr, 0)
        for j in range(len(transmissions)):
            planet_matrix = updateTransmittance(planet_matrix, transmissions[j][1])
            transmissions[j][2][i] = np.sum(planet_matrix*star_matrix)
        plt.figure()
        plt.imshow(planet_matrix*star_matrix)
        plt.title(str(i+1))
    return transmissions


 
#This creates a mask for the star by assigning all values in the matrix where the condition that their distance from the center is less than equal to the radius is true.
def create_star_mask(matrix_input, center, radius, mu):  
    if mu < 0 or mu>=1: raise ValueError("Limb darkening coefficient must be between 0 and +1.")
    Y,X = np.ogrid[:matrix_input.shape[0],:matrix_input.shape[1]]
    cart_dist = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    cos_matrix = np.zeros(cart_dist.shape)
    cos_matrix[cart_dist <= radius] = 1-mu*((1-np.sqrt((radius**2 - cart_dist[cart_dist <= radius]**2)/(radius**2))))
    
    plt.figure()
    plt.imshow(cos_matrix)
    return cos_matrix
    



#This creates a mask for the planet with the transmittance of red
def create_planetary_mask(planet_matrix, center, planetary_radius, annular_radius, transmittance):
    if annular_radius < 0: raise ValueError("Must have a non-negative atmospheric radius")
    elif annular_radius == 0: transmittance = 0
    elif transmittance <= 0.00 or transmittance >=1.00: raise ValueError("Transmittance must be between 0 and 1")
    if(annular_radius == 0 and transmittance!=0):
        print("Are you sure you provided your simulated exoplanet with an atmosphere?")
    
    Y, X = np.ogrid[:planet_matrix.shape[0], :planet_matrix.shape[1]]
    dist = (np.sqrt((X - center[0])**2 + (Y-center[1])**2))
        
    #-- Make a dedicated mask that indicates pixels that are in the atmosphere anulus
    planet_matrix[dist<=(annular_radius+planetary_radius)] = transmittance
    planet_matrix[dist<=planetary_radius] = 0
    
    
    
    return (planet_matrix)




#For the transit simulator, make it so that everything works. Get the atmospheric mask to work and same for updating transmittance
#Then, for the N-body, put the transit parameters into the particle class.
#implement some form of distance scaling and use the x-y coord stuff and get that working.   


def runsim(bodies, nsteps, tstep, transitActivated):
    XPosArr = np.full((len(bodies),nsteps), np.nan)
    YPosArr = np.full((len(bodies),nsteps), np.nan)
    
           
            
    N = 2500
    stellar_radius = 650
    
    if(stellar_radius>=0.75*N): raise ValueError("Must have smaller stellar radius proportional to resolution")
    mu = 0
    star_matrix= create_star_mask(np.zeros((N,N)), [N/2 + 1, N/2 + 1], stellar_radius, mu)
            
    planetary_radius = 30
    annular_radius = 2
                
    transmissions = (["red", float(input("Give a red transmission value:")), np.full(25, np.nan)],
                                 ["green", float(input("Give a green transmission value:")), np.full(25, np.nan)],
                                 ["blue", float(input("Give a blue transmission value:")), np.full(25, np.nan)],
                                 ["grey", float(input("Give a visual transmission value:")), np.full(25, np.nan)],
                                 ["yellow", float(input("Give a ultraviolet transmission value:")), np.full(25, np.nan)],
                                 ["orange", float(input("Give a infrared transmission value:")), np.full(25, np.nan)])
    for steps in range(nsteps):
        for i in range(len(bodies)):
            for j in range(len(bodies)):
                if i!=j:
                    bodies[i].getForce(bodies[j])
                    bodies[j].getForce(bodies[i])
            bodies[i].update(tstep)
        for l in range(len(bodies)):
            bodies[l].resetForce()
            XPosArr[l][steps] = bodies[l].x
            YPosArr[l][steps] = bodies[l].y
            if(XPosArr[l][steps]>0 and (YPosArr[l][steps]>0 and YPosArr[l][steps-1]<0)):
                if(transitActivated): 
                    hShift = -900
                    
                    planet_matrix = create_planetary_mask(np.ones((N, N)), [N/2 + hShift, N/2], planetary_radius, annular_radius, transmissions[0][1])

                    IntegratedSim  = runTransit(hShift, 72, 0, 0, 25, star_matrix, planet_matrix, transmissions)
                
            
                    plt.figure()
                    for b in range(len(transmissions)):
                        plt.plot(np.arange(25), normalize(IntegratedSim[b][2], star_matrix), color = IntegratedSim[b][0])
                    transitActivated = False
    return (XPosArr, YPosArr)

def plotter(XPosArr, YPosArr):
    plt.figure(figsize=(85,85))
    
    for ind in range(XPosArr.shape[0]):
        if ind == 0:
            plt.plot(XPosArr[ind,:]/150e9, YPosArr[ind,:]/150e9, "o", label='star', color = "orange")
        else:
            plt.plot(XPosArr[ind,:]/150e9, YPosArr[ind,:]/150e9, "-o", label='planet %d'%ind)#, color = "blue")
def setUpNBody(arr):
    #(StellarMass, PlanetMass, PlanetarySemiMajorAxis, PlanetVelo)
    Star = Particle(mass = arr[0], x = 0, y = 0, vx = 0, vy = 0)
    Planet = Particle(mass = arr[1], x = arr[2], y = 0, vx = 0, vy = arr[3])
    
    return [Star, Planet]

def run(tstep):
    plt.ion()
    params = inputPlanet()
    
    bodies = setUpNBody(params)
    
    steps = int(params[-1]/tstep)
    
    print("Steps = " + str(steps))
    
    XPosArr, YPosArr = runsim(bodies, tstep, steps,
                              True)
    plotter(XPosArr, YPosArr)
    
    
    
