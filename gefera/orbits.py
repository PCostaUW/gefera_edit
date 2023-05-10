import numpy as np
from timeit import default_timer as timer
import ctypes
import os
import fnmatch

__all__ = ['PrimaryOrbit', 'SatelliteOrbit', 'ConfocalOrbit']

class Orbit:
    
    """Parent class of all orbits. This class does not contain
    all the necessary attributes for all orbits. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    
    def __init__(self, a, t, e, p, w, i):
        
        self.a = a
        self.t = t
        self.e = e
        self.p = p
        self.w = w
        self.i = i
        
    def pdict(self):
        
        return vars(self)
    
class PrimaryOrbit(Orbit):
    
    """
    A heliocentric orbit for the primary body in the system. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    def def_e(ecosw,esinw):
        return np.sqrt(ecosw**2+esinw**2)
    def def_a(p1,b1,T,ecosw,esinw):
        return np.sqrt((1-b1**2)/(np.sin(Pi*T/p1*(1+esinw)**2/(1-def_e(ecosw,esinw)**2)**(3/2))**2)+b1**2)*(1+esinw)/(1-def_e(ecosw,esinw)**2)
    def def_i1(b1,p1,T,ecosw,esinw):
        return np.arccos(b1/def_a(p1,b1,T,ecosw,esinw)*(1+esinw)/(1-def_e(ecosw,esinw)**2))
    def def_w(ecosw,esinw):
        return np.arctan(esinw/ecosw)
    def def_t1(t0,p,ecosw,esinw):
        return t0+p*np.sqrt(1-def_e(ecosw,esinw)**2)/(2*Pi)*(def_e(ecosw,esinw)*np.sin(1.5*Pi-def_w(ecosw,esinw))/(1+def_e(ecosw,esinw)*np.cos(1.5*Pi-def_w(ecosw,esinw))) - 2/np.sqrt(1-def_e(ecosw,esinw)**2)*math.atan2(np.sqrt(1-def_e(ecosw,esinw)**2)*np.tan(0.75*Pi-0.5*def_w(ecosw,esinw)),(1+def_e(ecosw,esinw))))
    
    def __init__(self, T, t0, esinw, ecosw, p, b):
        a = def_a(p,b,T,ecosw,esinw)
        t = def_t1(t0,p,ecosw,esinw)
        e = def_e(ecosw,esinw)
        w = def_w(ecosw,esinw)
        super().__init__(a, t, e, p, w, i)
        
    def pdict(self):
        
        return {k + '1': v for k, v in vars(self).items()}
    
class SatelliteOrbit(Orbit):
    
    """
    The orbit of the moon around the planet.
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        o: Longitude of ascending node (in radians)
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
        m: Moon/planet mass ratio
    """
    #definitions of old parameter set in terms of new parameters
    def def_e(ecosw,esinw):
        return np.sqrt(ecosw**2+esinw**2)
    def def_t2(phi,p2):
        return phi*p2/(2*Pi)
    def def_w(ecosw,esinw):
        return np.arctan(esinw/ecosw)
    def def_i2(b2,a2,ecosw,esinw):
        return np.arccos(b2/a2*(1+esinw)/(1-def_e(ecosw,esinw)**2))
    
    def __init__(self, a, phi, esinw, ecosw, o, p, b, m):
        e=def_e(ecosw,esinw)
        t=def_t2(phi,p)
        w=def_w(ecosw,esinw)
        i=def_i2(b,a,ecosw,esinw)
        super().__init__(a, t, e, p, w, i)
        self.o = o
        self.m = m
        
    def pdict(self):
        
        return {k + '2': v for k, v in vars(self).items()}
        
class ConfocalOrbit(Orbit):
    
    """
    A second heliocentric orbit. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        o: Longitude of ascending node (in radians)
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    
    def __init__(self, a, t, e, p, o, w, i):
        super().__init__(a, t, e, p, w, i)
        self.o = o
        
    def pdict(self):
        
        return {k + '2': v for k, v in vars(self).items()}