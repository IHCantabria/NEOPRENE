'''
Library containing the mathematical functions of the Neyman-Scott point model Rectangular Pulses (NSRP)

	Authors: 
        + Javier Díez Sierra
	    + Salvador Navas Fernández
        + Manuel del Jesus
        
        
Point Model (NSRP)

	The NSRP process has the following parameters with their respective definitions:

	    landa --> Storms occur in a Poisson porcess of rate landa. Unidades: (time-1)
		mu_c  --> C number of cell origins is associated with each storm origin. An 
			  approach that ensures at least one rain cell follows each storm origin, 
			  is to take C-1 as a Poisson random variable with mean mu_c -1. 
			  Unidades: (numero de celldas por tormenta)
		eta   --> The duration of the pulse (rain cell) is an independent exponencial 
			  random variable with parameter eta. Unidades (time-1)
		beta  --> The cell origins being independently displaced from the storm origins
			  by distances that are exponetially distributed with parameter betha.
			  Unidades (time-1)
		xi    --> The intensity of the pulse can be exponencial, gamma or weibull 
			  distributed with parameteres escale pareameter xi. Unidades (mm-1)
		alpha --> the intensity of the pulse can be exponencial, gamma or weibull 
			  distributed with parameteres shape pareameter alpha. No tiene unidades.
'''


import numpy as np
import scipy as sp
import math as mt
from math import *


def E_X(xi, alpha, r):
    """Returns the intensity function (mm) that I am going to use (cell intensity).
       E' (exponential), 'W' (Weibull) and 'G' (Gamma).

    Exponential (1 parameters xi). The exponential function 
    is the same as the Weibull only with alpha=1. 
    Weibull (2 parameters, xi (scale parameter) and alpha (shape parameter).
    Gamma (2 parameters, xi (scale parameter) and alpha (shape parameter)
    
    xi:    Scale parameter
    alpha: Shape parameter
    r:     Exponent to which it is raised (e.g. E(X^r))"""
    
    E_X=(xi**(1*r))*mt.gamma(1+r*1)
    return E_X



def NSRP_mean(h, landa, mu_c, eta, xi, alpha, alpha_p):
    """Returns the average precipitation (mm/h) at a point.
    h : selected time interval (hours or days)"""

    result=list()
    for i in range(np.size(alpha_p)):
        result.append(h*landa[i]*mu_c[i]/(eta[i]*E_X(xi[i], alpha[i], 1)))
    return np.sum(result)

def NSRP_mean_ST(h, landa, mu_c, eta, xi, alpha, alpha_p):
    """Returns the average precipitation (mm/h) at a point.
    h : selected time interval (hours or days)"""

    result=list()
    for i in range(np.size(alpha_p)):
        result.append(h*landa[i]*mu_c[i]*E_X(xi[i], alpha[i], 1)/(eta[i]))
    return np.sum(result)


def NSRP_covariance(h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p):
    """Returns the covariance.
    Eq (12)--> Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0:
            A=eta[i]*h-1+exp(-eta[i]*h)
            B=betha[i]*h-1+exp(-betha[i]*h)
        elif l>0:
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1))
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1))
        result.append(alpha_p[i]*(landa[i]*(eta[i]**-3)*A*\
        (2*mu_c[i]*E_X(1/xi[i], alpha[i], 2) +\
	(E_X( 1/xi[i], alpha[i], 1)**2)*(betha [i]**2)*\
	(mu_c[i]**2)/((betha[i]**2)-(eta[i]**2)))-\
        landa[i]*(E_X(1/xi[i], alpha[i], 1)**2)*B*(mu_c[i]**2)/\
	(betha[i]*((betha[i]**2)-(eta[i]**2)))))
    return np.sum(result)
	

def multiply(numbers):  
    total = 1
    for x in numbers:
        total *= x  
    return total  


def NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p):
    """Returns the probability that it does not rain in an arbitrary time of length h.
    Eq (6) Stochastic point porcess modelling of rainfall. I. 
    sigle-site fitting and validation (Cowperwait, 1995)
    We use the approximation to avoid having to solve the integral--> Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(alpha_p)):

        result.append((exp(-landa[i]*h + (landa[i]*(betha[i]**-1)*((mu_c[i]-1)**-1))*(1-exp(1 -mu_c[i]+(mu_c[i]-1)*exp(-betha[i]*h)))-\
        (landa[i]/betha[i])*(0.577 + np.log(((eta[i]/(eta[i]-betha[i]))-exp(-betha[i]*h))*(mu_c[i]))))))
    return multiply(result)

def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
    
    
def fi_h(Datos, h, temporal_resolution):
    if temporal_resolution == 'D':
        tt = 0.1
    elif temporal_resolution == 'h':    
        tt = 0.001
    return np.sum(moving_average(Datos.values, h)<=tt)/(np.sum(Datos.values>=0))


def NSRP_fi_DD(h, landa, ipsilon, eta, betha, alpha_p):
    """Returns the probability that it will not rain in an arbitrary time of 
    length h conditional on the following time h not raining either.
    Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    result=NSRP_pdry(2*h, landa, ipsilon, eta, betha, alpha_p)/\
        NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p)
    return result

def NSRP_fi_WW(h, landa, ipsilon, eta, betha, alpha_p):
    """Returns the probability that it will rain in an arbitrary time of 
    length h conditional on the following time h also being rained on."""
    result=(1-2*NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p) +\
                       NSRP_pdry(2*h, landa, ipsilon, eta, betha, alpha_p))/\
                       (1-NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p))
    return result


def fi_DD(Datos, h, temporal_resolution):
    """(7, 8, 9). Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    return fi_h(Datos, 2*h, temporal_resolution)/fi_h(Datos, h, temporal_resolution)
    
    
def fi_WW(Datos, h, temporal_resolution):
    """(7, 8, 9). Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    return (1 - 2*fi_h(Datos, h, temporal_resolution)+fi_h(Datos, 2*h, temporal_resolution))/(1-fi_h(Datos, h, temporal_resolution))


def f_eta_beta_h(neta, beta, h):
    f=-2*neta**3*beta**2*exp(-neta*h)-2*neta**3*beta**2*exp(-beta*h)+neta**2*beta**3*exp(-2*neta*h)+2*neta**4*beta*exp(-neta*h)+2*neta**4*beta*exp(-beta*h)\
    +2*neta**3*beta**2*exp(-(neta+beta)*h)-2*neta**4*beta*exp(-(neta+beta)*h)-8*neta**3*beta**3*h+11*neta**2*beta**3-2*neta**4*beta\
    +2*neta**3*beta**2+4*neta*beta**5*h+4*neta**5*beta*h-7*beta**5-4*neta**5+8*beta**5*exp(-neta*h)-beta**5*exp(-2*neta*h)\
    -2*h*neta**3*beta**3*exp(-neta*h)-12*neta**2*beta**3*exp(-neta*h)+2*h*neta*beta**5*exp(-neta*h)+4*neta**5*exp(-beta*h)
    return f
def g_eta_beta_h(neta, beta, h):
    g=12*neta**5*beta*exp(-beta*h)+9*neta**4*beta**2+12*neta*beta**5*exp(-neta*h)+9*neta**2*beta**4+12*neta**3*beta**3*exp(-(neta+beta)*h)\
    -neta**2*beta**4*exp(-2*neta*h)-12*neta**3*beta**3*exp(-beta*h)-9*neta**5*beta-9*neta*beta**5-3*neta*beta**5*exp(-2*neta*h)\
    -neta**4*beta**2*exp(-2*beta*h)-12*neta**3*beta**3*exp(-neta*h)+6*neta**5*beta**2*h-10*beta**4*neta**3*h+6*beta**5*neta**2*h\
    -10*beta**3*neta**4*h+4*beta**6*neta*h-8*beta**2*neta**4*exp(-beta*h)+4*beta*neta**6*h+12*beta**3*neta**3\
    -8*beta**4*neta**2*exp(-neta*h)-6*neta**6-6*beta**6-2*neta**6*exp(-2*beta*h)-2*beta**6*exp(-2*neta*h)\
    +8*neta**6*exp(-beta*h)+8*beta**6*exp(-neta*h)-3*beta*neta**5*exp(-2*beta*h)
    return g

def NSRP_moments_order_3(C,h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p):
    """Returns the third-order moments.
    Eq (12)--> Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0:
            A=eta[i]*h-1+exp(-eta[i]*h)
            B=betha[i]*h-1+exp(-betha[i]*h)
        elif l>0:
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1))
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1))
	
        mu_x = E_X( 1/xi[i], alpha[i], 1)
        Ex2  = E_X( 1/xi[i], alpha[i], 2) 
        Ex3  = E_X( 1/xi[i], alpha[i], 3)
        
        if C=='Poisson':
            E_C1=mu_c[i]**2; E_C2=mu_c[i]**3
        if C=='Geometric':
            E_C1=2*mu_c[i]*(mu_c[i]-1); E_C2=6*mu_c[i]*(mu_c[i]-1)**2
        
        result.append(((6*landa[i]*mu_c[i]*Ex3*(eta[i]*h-2+eta[i]*h*exp(-eta[i]*h)\
		+2*exp(-eta[i]*h))/(eta[i]**4))+(3*landa[i]*mu_x*Ex2*(E_C1)*\
		f_eta_beta_h(eta[i], betha[i],h))/(2*(eta[i]**4)*betha[i]*\
		((betha[i]**2)-(eta[i]**2))**2))+landa[i]*(mu_x**3)*(E_C2)*\
		g_eta_beta_h(eta[i], betha[i],h)/(2*(eta[i]**4)*betha[i]*((eta[i]**2)-\
		(betha[i]**2))*(eta[i]-betha[i])*(2*betha[i]+eta[i])*(betha[i]+2*eta[i])))
    return np.sum(result)


