#This file has functions that calculate the old parameter set from the new parameter set
#These functions are called after the user inputs the new parameter values when initializing an orbit using orbits.py
#These are needed because the orbit coordinates and flux are calculated using the old parameter set
import numpy as np

Pi = np.pi

def ecc(ecosw,esinw):
    '''Calculates the eccentricity of the orbit.
    
    Args:
        ecosw: eccentricity multiplied by cos of the argument of periastron w
        esinw: eccentricity multiplied by sin of the argument of periastron w'''
    return np.sqrt(ecosw**2+esinw**2)

def t2(phi,p2):
    '''Calculates the time of periastron passage of the moon about the planet.
    
    Args:
        phi: 2*Pi*t2/p2
        p2: The orbital period of the moon about the planet'''
    return phi*p2/(2*Pi)

def a1(T,b1,p1,ecosw,esinw):
    '''Calculates the semi-major axis of the planet's orbit.
    Args:
        T: transit duration
        b1: impact parameter of the planet
        p1: orbital period of the planet
        ecosw: eccentricity multiplied by cos of the argument of periastron w
        esinw: eccentricity multiplied by sin of the argument of periastron w
    '''
    return np.sqrt((1-b1**2)/(np.sin(Pi*T/p1*(1+esinw)**2/(1-ecc(ecosw,esinw)**2)**(3/2))**2)+b1**2)*(1+esinw)/(1-ecc(ecosw,esinw)**2)

def i1(b1,p1,T,ecosw,esinw):
    '''Calculates the inclination of the planet's orbit.
    Args:
        b1: impact parameter
        p1: orbital period
        T: transit duration
        ecosw: component of eccentricity
        esinw: component of eccentricity
    '''
    return np.arccos(b1/a1(p1,b1,T,ecosw,esinw)*(1+esinw)/(1-ecc(ecosw,esinw)**2))

def i2(b2,a2,ecosw,esinw):
    '''Calculates the inclination of the moon's orbit.
    Args:
        b2: impact parameter
        a2: semi-major axis
        ecosw: component of eccentricity
        esinw: component of eccentricity
    '''
    return np.arccos(b2/a2*(1+esinw)/(1-ecc(ecosw,esinw)**2))

def omega(ecosw,esinw):
    '''Calculates the argument of periastron from the compenents of eccentricity'''
    return np.arctan(esinw/ecosw)

def t1(t0,p,ecosw,esinw):
    '''Calculates the time of periastron passage of the planet.
    Args:
        t0: time of mid-transit
        p: orbital period of the planet
        ecosw: component of eccentricity
        esinw: component of eccentricity
    '''
    return t0+p*np.sqrt(1-ecc(ecosw,esinw)**2)/(2*Pi)*(ecc(ecosw,esinw)*np.sin(1.5*Pi-omega(ecosw,esinw))/(1+ecc(ecosw,esinw)*np.cos(1.5*Pi-omega(ecosw,esinw))) - 2/np.sqrt(1-ecc(ecosw,esinw)**2)*math.atan2(np.sqrt(1-ecc(ecosw,esinw)**2)*np.tan(0.75*Pi-0.5*omega(ecosw,esinw)),(1+ecc(ecosw,esinw))))

