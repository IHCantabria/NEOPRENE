'''
Librería que contine las funciones matématicas del modelo puntual Neyman-Scott 
Rectangular Pulses (NSRP) y del modelo espacio-temporal Neyman-Scott Rectangular Pulses (STNSRP)

	Autor: Javier Díez Sierra
	Fechas 06-02-2017

Los estadísticos puntuales del modelo STNSRP son ajustados con las mismas ecuaciones con las que ajustamos el modelo puntual NSRP. 
Los estadísticos y parámetros propios del modelo espacial se comentaran más adelante.

a)Modelo Puntual (NSRP)

	El proceso de NSRP tiene los siguentes parámetros con sus respectivas definiciones:

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

	Existen dos variantes en el proceso de NSTP: 

		a) El modelo se ajusta para disponer de más de un celda tipo (el proceso no tiene
                   una única celda tipo) si no varias y cada una de ellas tiene diferentes 
                   parámetros de eta, xi y alpha.
		   Además, cada celda tipo tine una probabilidad de ocurrencia alpha_p. En este 
		   caso solamente esta programando para dos celdas tipo. La primera con una 
		   probabilidad de ocurrencia alpha_p y la segunda 1-alpha_p.(LA OPCIÓN PARA
		   SIMULAR EL MODELO NSRP PARA DIFERENTES CELDAS TIPO NO FUNCIONA YA QUE LAS 	
		   ECUACIONES SON DIFERENTES)

		b) Aquella que se basa en varios processos de Neyman-Scott superpuestos. En este   
                   caso se ajusta el modelo para que cada tormenta sea totalmente diferente luego   
                   todos sus parámetros serán independientes y alpha_p no se utila para nada.
		   (EN ESTE CASO COMO LOS PROCESO SON INDPENDIENTES LAS ECUACIONES SON LAS MISMAS Y 			   FUNCIONA PERFECTAMENTE)'''




'''NSRP--> Propiedades matemáticas del modelo puntual'''

import numpy as np
import scipy as sp
import math as mt
from math import *


def E_X(I_F, xi, alpha, r):
    """Devuelve funcion de intenidad (mm) que voy a utilizar (cell intensity)
    I_F: Funcion de intensidad que queremos ajustar. Podemos elegir entre: 
    'E' (exponencial), 'W' (Weibull) y 'G' (Gamma).(SOLO FUNCIONA PARA 'E')

    Exponencial (1 parámetros xi). La función exponencial 
    es igual que la Weibull solo que con alpha=1. 
    Weibull (2 parámetros, xi (parámetro escala) y alpha (parámetro forma)
    Gamma (2 parámetros, xi (parámetro escala) y alpha (parámetro forma)
    
    xi: Parámetro de escala
    alpha: Parámetro de forma
    r:Exponente al que está elevado (Ej. E(X^r))"""
    
    if I_F=='E': ##Exponential
        #E_X=(xi**(alpha*r))*mt.gamma(1+r*alpha)
        E_X=(xi**(1*r))*mt.gamma(1+r*1)
    elif I_F=='W': ##Weibull
        E_X=(xi**(alpha*r))*mt.gamma(1+r*alpha)
        #E_X=(xi**(r))*mt.gamma(1+r/alpha)#Mixed rectangular pulses models of rainfall
    elif I_F=='G': ##Gamma
        E_X=((xi**r)*mt.gamma(alpha + r))/mt.gamma(alpha)
    return E_X



def NSRP_mean(I_F, h, landa, mu_c, eta, xi, alpha, alpha_p):
    """Devuelve la precipitación media (mm/h) en un punto.
    h : intervalo de tiempo seleccionado (hours or days)"""

    result=list()
    for i in range(np.size(alpha_p)):
        result.append(h*landa[i]*mu_c[i]/(eta[i]*E_X(I_F, xi[i], alpha[i], 1)))
    return np.sum(result)

def NSRP_mean_ST(I_F, h, landa, mu_c, eta, xi, alpha, alpha_p):
    """Devuelve la precipitación media (mm/h) en un punto.
    h : intervalo de tiempo seleccionado (hours or days)"""

    result=list()
    for i in range(np.size(alpha_p)):
        result.append(h*landa[i]*mu_c[i]*E_X(I_F, xi[i], alpha[i], 1)/(eta[i]))
    return np.sum(result)




def NSRP_covariance(I_F, h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p):
    """Devuelve la covarianza.
    Eq (12)--> Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0:#Varianza
            A=eta[i]*h-1+exp(-eta[i]*h);#A_h_0
            B=betha[i]*h-1+exp(-betha[i]*h)#B_h_0
        elif l>0:#Temporal Covarianza
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1))#A_h_l
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1))#B_h_l
        #Eq A Further Developments of the Neyman-Scott Clustered Point 
        #Process for Modeling Rainfall  (Cowperwait, 1991)
        #result.append((landa[i]*((ipsilon[i]**2)-1)*((beta[i]**3)*
        #A-(eta[i]**3)*B)/(beta[i]*(xi [i]**2)*\
        #(eta[i]**3)*((beta[i]**2)-(eta[i]**2))))+4*landa[i]*ipsilon[i]*A/((xi[i]**2)*(eta[i]**3)))
        #Eq B
        #Ecuacion: A space-time Neyman-Scott model of rainfall: 
        #Empirical anlysis of extrems (Cowperwait, 2002)
        result.append(alpha_p[i]*(landa[i]*(eta[i]**-3)*A*\
        (2*mu_c[i]*E_X(I_F, 1/xi[i], alpha[i], 2) +\
	(E_X(I_F, 1/xi[i], alpha[i], 1)**2)*(betha [i]**2)*\
	(mu_c[i]**2)/((betha[i]**2)-(eta[i]**2)))-\
        landa[i]*(E_X(I_F, 1/xi[i], alpha[i], 1)**2)*B*(mu_c[i]**2)/\
	(betha[i]*((betha[i]**2)-(eta[i]**2)))))
    return np.sum(result)
	#Ecuacion: Some models for rainfall based based on stochastic point processes
	#result.append(landa[i]*(eta[i]**-3)*A*\
	#(2*ipsilon[i]*E_X(I_F, xi[i], alpha[i], 2) + (E_X(I_F, xi[i], alpha[i], 1)**2)*(beta[i]**2)*((ipsilon[i])*(ipsilon[i]+2))/((beta[i]**2)-(eta[i]**2)))-\
	#landa[i]*(E_X(I_F, xi[i], alpha[i], 1)**2)*B*((ipsilon[i])*(ipsilon[i]+2))/(beta[i]*((beta[i]**2)-(eta[i]**2))))
	#Ecuacuines del puntual son iguales a las de el espacial
	

def multiply(numbers):  
    total = 1
    for x in numbers:
        total *= x  
    return total  


def NSRP_pdry(h, landa, mu_c, eta, betha, alpha_p):
    """Devuelve la probabilidad de que no llueva en un tiempo arbitrario de longitud h.
    Eq (6) Stochastic point porcess modelling of rainfall. I. 
    sigle-site fitting and validation (Cowperwait, 1995)
    Usamos la aproximación para no tener que resolver la intregral--> 
    Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    ##Eq (6) Stochastic point porcess modelling of rainfall. I. sigle-site fitting and validation (Cowperwait, 1995)
    ## Aproximación intregral--> Further Developments of the Neyman-Scott Clustered Point Process for Modeling Rainfall (Cowperwait, 1991)
    result=list()
    for i in range(np.size(alpha_p)):
        #result.append(exp(-landa[i]*h + (landa[i]*(betha[i]**-1)*((mu_c[i])**-1))*(1-exp(-mu_c[i]+(mu_c[i])*exp(-betha[i]*h)))-\
        #   (landa[i]/betha[i])*(0.577 + log(((eta[i]/(eta[i]-betha[i]))-exp(-betha[i]*h))*(mu_c[i]-1)))))##ojo he puesto el último -1

        result.append((exp(-landa[i]*h + (landa[i]*(betha[i]**-1)*((mu_c[i]-1)**-1))*(1-exp(1 -mu_c[i]+(mu_c[i]-1)*exp(-betha[i]*h)))-\
        (landa[i]/betha[i])*(0.577 + np.log(((eta[i]/(eta[i]-betha[i]))-exp(-betha[i]*h))*(mu_c[i]))))))##ojo he puesto el último -1 asi funciona con exponencial
    return multiply(result)##OJO The probability that an arbitrary time interval is dry at a point is obtained
    #by multiplying the probabilities of the independent porcesses (Basque Country, Spain. P.Cowpertwait)


#Veo que da diferente y me quedo coen esta 
def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
def fi_h(Datos, h):
    ##Calcula la probabilidad de que no llueva in un tiempo arbitrario de longitud h
    return np.sum(moving_average(Datos.values, h)<=0.1)/(np.sum(Datos.values>=0))
    ##MUY IMPORTANTE METER EK <=0.1 PORQUE EL MODELO SIMEPRE DA LLUVIA Y SI NO LO PONEMOS EL PDRY DARIA 0


def NSRP_fi_DD(h, landa, ipsilon, eta, betha, alpha_p):
    """Devuelve la probabilidad de que no llueva in un tiempo arbitrario de 
    longitud h condicionado a que el siguiente tiempo h tampoco llueva
    Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    result=NSRP_pdry(2*h, landa, ipsilon, eta, betha, alpha_p)/\
        NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p)
    return result

def NSRP_fi_WW(h, landa, ipsilon, eta, betha, alpha_p):
    """Devuelve la probabilidad de que llueva in un tiempo arbitrario de 
    longitud h condicionado a que el siguiente tiempo h también llueva."""
    result=(1-2*NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p) +\
                       NSRP_pdry(2*h, landa, ipsilon, eta, betha, alpha_p))/\
                       (1-NSRP_pdry(h, landa, ipsilon, eta, betha, alpha_p))
    return result


#def fi_DD(Datos, h):
#    """Calcula la probabilidad de que no llueva in un tiempo arbitrario de 
#    longitud h condicionado a  que el siguiente tiempo h tampoco llueva"""
#
#    n=0
#    for i in range(len(Datos)): 
#        if np.nansum(Datos.values[i-h:i+h])==0: n=n+1   
#    return n/(np.sum(moving_average(Datos.values, h)==0))


#def fi_WW(Datos, h):
#    """Calcula la probabilidad de que llueva in un tiempo arbitrario de longitud 
#    h condicionado a que el siguiente tiempo h llueva"""
#    n=0
#    for i in range(len(Datos)): 
#        if np.nansum(Datos.values[i-h:i]>0) and np.nansum(Datos.values[i:i+h])>0: n=n+1   
#    return n/(np.sum(moving_average(Datos.values, h)>0))

##Transition probabilities
def fi_DD(Datos, h):
    """(7, 8, 9). Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    return fi_h(Datos, 2*h)/fi_h(Datos, h)
def fi_WW(Datos, h):
    """(7, 8, 9). Cowperwait , 1995. Stochastic point process modelling of rainfall. I. 
    single-site fitting and validation. NSRP"""
    return (1 - 2*fi_h(Datos, h)+fi_h(Datos, 2*h))/(1-fi_h(Datos, h))

##Momentos de tercer orden
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

def NSRP_momentos_orden_3(I_F,C,h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p):
    """Devuelve los momentos de tercer orden.
    Eq (12)--> Further Developments of the Neyman-Scott Clustered Point 
    Process for Modeling Rainfall (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0:#Varianza
            A=eta[i]*h-1+exp(-eta[i]*h);#A_h_0
            B=betha[i]*h-1+exp(-betha[i]*h)#B_h_0
        elif l>0:#Temporal Covarianza
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1))#A_h_l
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1))#B_h_l  
	
        mu_x = E_X(I_F, 1/xi[i], alpha[i], 1)
        Ex2  = E_X(I_F, 1/xi[i], alpha[i], 2) 
        Ex3  = E_X(I_F, 1/xi[i], alpha[i], 3)
        
        if C=='Poisson':
            E_C1=mu_c[i]**2; E_C2=mu_c[i]**3
            #E_C1=mu_c[i]*(mu_c[i]+2); ##Iturbe: Some models..
        if C=='Geometric':
            E_C1=2*mu_c[i]*(mu_c[i]-1); E_C2=6*mu_c[i]*(mu_c[i]-1)**2
        
        result.append(((6*landa[i]*mu_c[i]*Ex3*(eta[i]*h-2+eta[i]*h*exp(-eta[i]*h)\
		+2*exp(-eta[i]*h))/(eta[i]**4))+(3*landa[i]*mu_x*Ex2*(E_C1)*\
		f_eta_beta_h(eta[i], betha[i],h))/(2*(eta[i]**4)*betha[i]*\
		((betha[i]**2)-(eta[i]**2))**2))+landa[i]*(mu_x**3)*(E_C2)*\
		g_eta_beta_h(eta[i], betha[i],h)/(2*(eta[i]**4)*betha[i]*((eta[i]**2)-\
		(betha[i]**2))*(eta[i]-betha[i])*(2*betha[i]+eta[i])*(betha[i]+2*eta[i])))
    return np.sum(result)

'''
b) STNSRP

	En el modelo espacial además de las propiedades puntuales ajustamos las espaciales  
        (correlación espacial). Como mínimo aparece un parámero nuevo que es el radio de la celda pero
	también podemos incluir el radio de la tormenta que para grades cuencas mejora la 
        correlación espacial haciendo que no llueva en toda la cuenca a la vez.
	
        Además de los parámetros definidos en el modelo NSRP, el modelo STNSRP tiene los siguentes    
        parámetros con sus respectivas definiciones:
        
        fi_may  --> the cell radius is a independent exponencial random variable with parameter  
                    fi_may. Unidades (km-1). 
        fi_may_s--> the storm radius is a independent exponencial random variable with parameter  
                    fi_may_s. Unidades (km-1).
        fi_min  --> En este caso, una vez calculado el parámetro ipsilon-1=mu_c, habría que calcular el  
                    parámetro fi_min spatial rate (number of cells per storm) is a independent 
                    exponencial random variable with parameter fi_min. Unidades (cells/km^2). fi_min  
                    depende exclusivamente de fi_may y ipsilon'''


'''STNSRP--> Propiedades matemáticas del modelo espacial'''

def STNSRP_fi_min(mu_c, fi_may):
    """Devuelve el numero de celdas por quilómetro cuadrado. 
    Ecuación 3: A space-time Neyman-Scott model of rainfall: Empirical analysis of extremes"""

    return (mu_c*(fi_may**2))/(2*np.pi)

def NSRP_cross_correlation(I_F, h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p,fi_may, d):
    """Devuleve la correlación cruzada
    Eq (12)--> Further Developments of the Neyman-Scott Clustered Point Process 
    for Modeling Rainfall  (Cowperwait, 1991)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0: ##Varianza
            A=eta[i]*h-1+exp(-eta[i]*h); #A_h_0
            B=betha[i]*h-1+exp(-betha[i]*h) #B_h_0
        elif l>0: ##Temporal Covarianza
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1)) #A_h_l
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1)) #B_h_l

        result.append(NSRP_covariance(I_F, h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
        -2*landa[i]*(1-P_fi_d(fi_may[i], d))*mu_c[i]*E_X(I_F, 1/xi[i], alpha[i], 2)*A/(eta[i]**3))

    return multiply(result)


def NSRP_cross_correlation_with_storm_radius(I_F, h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p, fi_may, d, fi_may_storm):
    """Devuleve la correlación cruzada (incluye el radio de la tormenta)
    Eq (A3)--> Regionalised spatiotemporal rainfall and temperature models for 
    flood studies in the Basque Country, Spain (Cowpertwait, 2012)"""

    result=list()
    for i in range(np.size(landa)):
        if l==0: ##Varianza
            A=eta[i]*h-1+exp(-eta[i]*h); #A_h_0
            B=betha[i]*h-1+exp(-betha[i]*h) #B_h_0
        elif l>0: ##Temporal Covarianza
            A=0.5*((1-exp(-eta[i]*h))**2)*exp(-eta[i]*h*(l-1)) #A_h_l
            B=0.5*((1-exp(-betha[i]*h))**2)*exp(-betha[i]*h*(l-1)) #B_h_l

        result.append(P_fi_d(fi_may_storm[i], d)*NSRP_covariance(I_F, h,l, landa, mu_c, eta, xi, betha, alpha, alpha_p)\
        -2*landa[i]*(1-P_fi_d(fi_may[i], d))*mu_c[i]*E_X(I_F, 1/xi[i], alpha[i], 2)*A/(eta[i]**3))

    return np.sum(result)

def fy(fi_may, d, y):
    result=((fi_may*d/(2*cos(y)))+1)*exp(-fi_may*d/(2*cos(y)))
    return result

def sum_fy(fi_may, d):
    sum=0
    for ii, i in enumerate(np.arange(1,5)):
        aux=2*fy(fi_may, d, 2*np.pi*i/20) + 4*fy(fi_may, d, (((2*np.pi*i)+np.pi)/20))
        sum=sum+aux
    return sum
        
def P_fi_d(fi_may, d):
    """Probability that a cell overlaps a point x given that it overlapped a point y a distance d from x.     
    Ecuación 9: A space-time Neyman-Scott model of rainfall: Empirical analysis of extremes"""

    result=(1/30)*sum_fy(fi_may, d) -(1/30)*fy(fi_may, d, 0.0)
    return result




