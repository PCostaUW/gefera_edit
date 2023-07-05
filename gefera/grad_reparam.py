import numpy as np
from .reparam import *

Pi = np.pi

#new derivatives of T
def dTdp(P,a,w,e,i):
    return 1/Pi * (1-e**2)**1.5 / (1+e*np.sin(w))**2 *np.arcsin(np.sqrt(1-a**2*(1-e**2)**2/(1+e*np.sin(w))**2*np.cos(i)**2)/(a*(1-e**2)/(1+e*np.sin(w))*np.sin(i)))

def dTda(P,a,w,e,i):
    b = P/Pi * (1-e**2)**(3/2)/(1+e*np.sin(w))**2
    c = (1-e**2)**2/(1+e*np.sin(w))**2*np.cos(i)**2
    d = (1-e**2)/(1+e*np.sin(w))*np.sin(i)
    
    return -b/(d*a**2 *np.sqrt(1 - c* a**2) *np.sqrt(c/d**2 - 1/(d**2 * a**2) + 1))

def dTdw(P,a,w,e,i):
    b = P/Pi*(1-e**2)**(3/2)
    c = a**2*(1-e**2)**2*np.cos(i)**2
    d = a*(1-e**2)*np.sin(i)
    return (e*b*np.cos(w)*((1 + e *np.sin(w))/(d*np.sqrt(1 - c/(1 + e*np.sin(w))**2)*np.sqrt((c + d**2 - e**2 *np.sin(w)**2 - 2*e*np.sin(w) - 1)/d**2)) - 2*np.arcsin(((1 + e*np.sin(w)) *np.sqrt(1 - c/(1 + e *np.sin(w))**2))/d)))/(1 + e*np.sin(w))**3

def dTde(P,a,w,e,i):
    q = P/Pi
    b = a**2*np.cos(i)**2
    d = a*np.sin(i)
    return -((2 *(1-e**2)**(3/2)* q* np.arcsin(((1+e *np.sin(w))* np.sqrt(1-(b *(1-e**2)**2)/(1+e*np.sin(w))**2))/(d *(1-e**2)))* np.sin(w))/(1+e*np.sin(w))**3)-(3*e*np.sqrt(1-e**2)* q *np.arcsin(((1+e *np.sin(w)) *np.sqrt(1-(b* (1-e**2)**2)/(1+e *np.sin(w))**2))/(d *(1-e**2))))/(1+e *np.sin(w))**2+((1-e**2)**(3/2)* q *(((1+e *np.sin(w))* ((2 *b *(1-e**2)**2 *np.sin(w))/(1+e*np.sin(w))**3+(4 *b *e *(1-e**2))/(1+e *np.sin(w))**2))/(2*d*(1-e**2)* np.sqrt(1-(b *(1-e**2)**2)/(1+e* np.sin(w))**2))+(np.sin(w)* np.sqrt(1-(b *(1-e**2)**2)/(1+e *np.sin(w))**2))/(d *(1-e**2))+(2 *e *(1+e *np.sin(w))* np.sqrt(1-(b *(1-e**2)**2)/(1+e *np.sin(w))**2))/(d*(1-e**2)**2)))/((1+e* np.sin(w))**2 *np.sqrt(1-((1+e *np.sin(w))**2 *(1-(b*(1-e**2)**2)/(1+e *np.sin(w))**2))/(d**2 *(1-e**2)**2)))

def dTdi(P,a,w,e,i):
    b = P/Pi * (1-e**2)**(3/2)/(1+e*np.sin(w))**2
    c = a**2 * (1-e**2)**2 / (1+e*np.sin(w))**2
    d = a * (1-e**2)/(1+e*np.sin(w))
    return (b *((c *np.cos(i))/(d *np.sqrt(1-c*np.cos(i)**2))-(np.sqrt(1-c *np.cos(i)**2) *np.cos(i))/(np.sin(i)**2*d)))/np.sqrt(1-((1-c *np.cos(i)**2))/(np.sin(i)**2*d**2))

#derivatives of t0
def dt0dt1():
    return 1
def dt0de(t,e,p,w):
    return (e *p *(-((2 *np.arctan((np.sqrt(1-e**2) *np.tan(1/2 *(1.5*Pi -w)))/(1+e)))/np.sqrt(1-e**2))+(e *np.sin(1.5*Pi -w))/(1+e *np.cos(1.5*Pi -w))))/(2 *np.sqrt(1-e**2)* Pi)-(np.sqrt(1-e**2)* p *(-((2 *e *np.arctan((np.sqrt(1-e**2) *np.tan(1/2 *(1.5*Pi -w)))/(1+e)))/(1-e**2)**(3/2))-(e*np.cos(1.5*Pi -w) *np.sin(1.5*Pi -w))/(1+e *np.cos(1.5*Pi -w))**2+np.sin(1.5*Pi -w)/(1+e *np.cos(1.5*Pi -w))-(2 *(-((e *np.tan(1/2 *(1.5*Pi -w)))/((1+e) *np.sqrt(1-e**2)))-(np.sqrt(1-e**2) *np.tan(1/2 *(1.5*Pi -w)))/(1+e)**2))/(np.sqrt(1-e**2)* (1+((1-e**2) *np.tan(1/2 *(1.5*Pi -w))**2)/(1+e)**2))))/(2 *Pi)
def dt0dp(t,e,p,w):
    return -np.sqrt(1-e**2) /2/Pi * (e*np.sin(1.5*Pi-w)/(1+e*np.cos(1.5*Pi-w)) - 2/np.sqrt(1-e**2)*np.arctan(np.sqrt(1-e**2)*np.tan(0.75*Pi-0.5*w)/(1+e)))
def dt0dw(t,e,p,w):
    return -((np.sqrt(1-e**2) *p* (-((e *np.cos(1.5*Pi -w))/(1+e *np.cos(1.5*Pi -w)))-(e**2 *np.sin(1.5*Pi -w)**2)/(1+e *np.cos(1.5*Pi -w))**2+1/np.cos(1/2 *(1.5*Pi -w))**2/((1+e) *(1+((1-e**2) *np.tan(1/2* (1.5*Pi -w))**2)/(1+e)**2))))/(2*Pi))

#derivatives of b
def dbda(a,i,e,w):
    return np.cos(i)*(1-e**2)/(1+e*np.sin(w))
def dbdi(a,i,e,w):
    return -a*np.sin(i)*(1-e**2)/(1+e*np.sin(w))
def dbde(a,i,e,w):
    return (-a*np.cos(i)*((e**2+1)*np.sin(w) + 2*e))/(1+e*np.sin(w))**2
def dbdw(a,i,e,w):
    return (a*np.cos(i))*e*(e**2-1)*np.cos(w)/(e*np.sin(w)+1)**2

#derivatives of phi
def dphidt(t,e,p,w):
    return 2*Pi/p
def dphide(t,e,p,w):
    return 2*Pi/p*(2 *e *np.arctan(((np.sqrt(1-e**2)* np.tan(1/2 *(3*Pi/2 -w))))/(1+e))/(1-e**2)**(3/2)+np.sqrt(1-e**2) *p *np.cos(3*Pi/2 -w) *np.sin(3*Pi/2 -w)/(2 *(1+e *np.cos(3*Pi/2 -w))**2)+e *p *np.sin(3*Pi/2 -w)/(2 *np.sqrt(1-e**2) *(1+e *np.cos(3*Pi/2 -w)))-(p/2*np.sin(1.5*Pi-w)*np.sqrt(1-e**2))/(1+e*np.cos(1.5*Pi-w))+2*(-e*np.tan(0.75*Pi-w/2)/((np.sqrt(1-e**2))*(1+e))-np.sqrt(1-e**2)*np.tan(0.75*Pi-w/2)/(1+e)**2)/(np.sqrt(1-e**2)*(1+np.tan(0.75*Pi-w/2)**2*(1-e**2)/(1+e)**2))+2*np.arctan(np.sqrt(1-e**2)/(1+e)*np.tan(0.75*Pi-w/2))/np.sqrt(1-e**2))
def dphidp(t,e,p,w):
    esinw = e*np.sin(w)
    ecosw = e*np.cos(w)
    return (-1/p**2)*2*Pi*t0(e,w,p,t) + 2*Pi/p*dtpdp(p,esinw,ecosw)
def dphidw(t,e,p,w):
    return 2*Pi/p*((p*np.sqrt(1-e**2)*np.cos(3*Pi/2-w))/(2*(1+e*np.cos(3*Pi/2-w)))+(p*e*np.sqrt(1-e**2)*np.sin(3*Pi/2-w)**2)/((2*(1+e*np.cos(3*Pi/2-w)))**2)-1/(np.cos(3*Pi/4-w/2)**2*(1+e)*(1+(1-e**2)/((1+e)**2)*np.tan(3*Pi/4-w/2))))

#derivatives of sigma = ecosw and rho = esinw
def dsigde(e,w):
    return np.cos(w)
def dsigdw(e,w):
    return -e*np.sin(w)

def drhode(e,w):
    return np.sin(w)
def drhodw(e,w):
    return e*np.cos(w)

#derivative of old paramaters with respect to new
def di2db(a,esinw,ecosw,b):
    return -(1+esinw)/(a*(1-ecc(ecosw,esinw)**2)*np.sqrt(1-b**2*(1+esinw)**2/(1-ecc(ecosw,esinw)**2)**2/a**2))

def di2decosw(b,a,esinw,ecosw):
    c = ecosw
    s = esinw
    return -((2 *b *c *(1+s))/(a *(1-c**2-s**2)**2 *np.sqrt(1-(b**2 *(1+s)**2)/(a**2 *(1-c**2-s**2)**2))))

def di2desinw(b,a,esinw,ecosw):
    c = ecosw
    s = esinw
    return -(((2 *b *s *(1+s))/(a *(1-c**2-s**2)**2)+b/(a *(1-c**2-s**2)))/np.sqrt(1-(b**2 *(1+s)**2)/(a**2 *(1-c**2-s**2)**2)))

def dadesinw(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -(((1-b**2)* (1+s)* ((3 *Pi* s* (1+s)**2 *T)/(p *(1-c**2-s**2)**(5/2))+(2*Pi* (1+s)* T)/(p* (1-c**2-s**2)**(3/2)))* 1/np.tan((Pi* (1+s)**2 *T)/(p* (1-c**2-s**2)**(3/2))) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/((1-c**2-s**2) *np.sqrt(b**2+(1-b**2)* 1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)))+(2 *s *(1+s) *np.sqrt(b**2+(1-b**2)* 1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2))/(1-c**2-s**2)**2+np.sqrt(b**2+(1-b**2)* 1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/(1-c**2-s**2)

def dadecosw(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -((3 *(1-b**2) *c *Pi *(1+s)**3 *T *1/np.tan((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2))) *1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/(p *(1-c**2-s**2)**(7/2) *np.sqrt(b**2+(1-b**2) *1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)))+(2 *c *(1+s) *np.sqrt(b**2+(1-b**2) *1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2))/(1-c**2-s**2)**2

def dadp(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return ((1-b**2) *Pi* (1+s)**3 *T *1/np.tan((Pi* (1+s)**2 *T)/(p* (1-c**2-s**2)**(3/2))) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/(p**2 *(1-c**2-s**2)**(5/2) *np.sqrt(b**2+(1-b**2) *1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2))
    #return 2/3*p**(-1/3)
def di1desinw(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -((b *(1-b**2)* ((3 *Pi *s *(1+s)**2 *T)/(p *(1-c**2-s**2)**(5/2))+(2 *Pi* (1+s)* T)/(p *(1-c**2-s**2)**(3/2))) *1/np.tan((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2))) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/((b**2+(1-b**2) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)**(3/2) *np.sqrt(1-b**2/(b**2+(1-b**2) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2))))

def di1decosw(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -((3 *b *(1-b**2) *c *Pi *(1+s)**2 *T *1/np.tan((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))* 1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)/(p *(1-c**2-s**2)**(5/2) *(b**2+(1-b**2) *1/np.sin((Pi* (1+s)**2 *T)/(p* (1-c**2-s**2)**(3/2)))**2)**(3/2) *np.sqrt(1-b**2/(b**2+(1-b**2) *1/np.sin((Pi *(1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2))))

def di1dp(b,T,p,esinw,ecosw):
    s = esinw
    c = ecosw
    return (b *(1-b**2) *Pi* (1+s)**2 *T *1/np.tan((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2))) *1/np.sin((Pi* (1+s)**2 *T)/(p* (1-c**2-s**2)**(3/2)))**2)/(p**2 *(1-c**2-s**2)**(3/2)* (b**2+(1-b**2) *1/np.sin((Pi* (1+s)**2 *T)/(p *(1-c**2-s**2)**(3/2)))**2)**(3/2)* np.sqrt(1-b**2/(b**2+(1-b**2) *1/np.sin((Pi* (1+s)**2 *T)/(p* (1-c**2-s**2)**(3/2)))**2)))

def dtpdesinw(p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -((p *s *(-(c/(1-s))+(2* np.arctan((np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(1+np.sqrt(c**2+s**2))))/np.sqrt(1-c**2-s**2)))/(2 *Pi *np.sqrt(1-c**2-s**2)))+(p *np.sqrt(1-c**2-s**2) *(-(c/(1-s)**2)+(2 *s *np.arctan((np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(1+np.sqrt(c**2+s**2))))/(1-c**2-s**2)**(3/2)+(2 *((np.sqrt(1-c**2-s**2) *1/np.cos(Pi/4+1/2 *np.arctan(s/c))**2)/(2 *c *(1+s**2/c**2)* (1+np.sqrt(c**2+s**2)))-(s* np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(np.sqrt(c**2+s**2)* (1+np.sqrt(c**2+s**2))**2)-(s *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(np.sqrt(1-c**2-s**2) *(1+np.sqrt(c**2+s**2)))))/(np.sqrt(1-c**2-s**2) *(1+((1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c))**2)/(1+np.sqrt(c**2+s**2))**2))))/(2 *Pi)

def dtpdecosw(p,esinw,ecosw):
    s = esinw
    c = ecosw
    return -((c *p* (-(c/(1-s))+(2 *np.arctan((np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(1+np.sqrt(c**2+s**2))))/np.sqrt(1-c**2-s**2)))/(2* Pi *np.sqrt(1-c**2-s**2)))+(p *np.sqrt(1-c**2-s**2)* (-(1/(1-s))+(2 *c *np.arctan((np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(1+np.sqrt(c**2+s**2))))/(1-c**2-s**2)**(3/2)+(2 *(-((s* np.sqrt(1-c**2-s**2)* 1/np.cos(Pi/4+1/2 *np.arctan(s/c))**2)/(2* c**2 *(1+s**2/c**2)* (1+np.sqrt(c**2+s**2))))-(c* np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(np.sqrt(c**2+s**2) *(1+np.sqrt(c**2+s**2))**2)-(c *np.tan(Pi/4+1/2 *np.arctan(s/c)))/(np.sqrt(1-c**2-s**2) *(1+np.sqrt(c**2+s**2)))))/(np.sqrt(1-c**2-s**2) *(1+((1-c**2-s**2) *np.tan(Pi/4+1/2 *np.arctan(s/c))**2)/(1+np.sqrt(c**2+s**2))**2))))/(2 *Pi)

def dtpdp(p,esinw,ecosw):
    c = ecosw
    s = esinw
    return (np.sqrt(1-c**2-s**2) *(-(c/(1-s))+(2 *np.arctan((np.sqrt(1-c**2-s**2) *np.tan(Pi/4+1/2 *omega(ecosw,esinw)))/(1+np.sqrt(c**2+s**2))))/np.sqrt(1-c**2-s**2)))/(2* Pi)

def de1desinw(esinw,ecosw):
    if ecosw <= 1e-20:
        return 1
    else:
        return 1/np.sqrt(1+np.tan(omega(ecosw,esinw))**2)

def de1decosw(esinw,ecosw):
    if ecosw <= 1e-20:
        return 1
    else:
        return 1/np.sqrt(1+1/np.tan(omega(ecosw,esinw))**2)

def dw1desinw(esinw,ecosw):
    return ecosw/(ecosw**2+esinw**2)

def dw1decosw(esinw,ecosw):
    return -esinw/(ecosw**2+esinw**2)