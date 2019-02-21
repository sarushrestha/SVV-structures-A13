# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 10:21:04 2019

@author: Saru Shrestha

"""
import numpy as np 
from math import *

#parameters
aileron_span = 2.691
x_1 = 0.174 
x_2 = 1.051
x_3 = 2.512
x_a1 = x_2 - 0.15
x_a2 =  x_2 + 0.15
aileron_deflection= (25.0/180.0) * pi 
q = 1000 #N/m

#variables
moment_inertia_zz = 1000 #m^4



def deflection(x):
    if 0<=x<=x_1:
        v = -(1/(moment_inertia_zz*e_modulus))*((-q*x^4)/24)
    elif x_1 <= x <= x_2:
        v = -(1/(moment_inertia_zz*e_modulus))*(r_1*(x^3/6 - x_1)- (q*x^4)/24)
    elif x_2 <= x <= x_3:
        v = -(1/(moment_inertia_zz*e_modulus))*(r_2* ( x^3/ 6 -x_2) + r_1*(x^3/6 - x_1)- (q*x^4)/24)
    elif x_3 <= x<= aileron_span:
        v = -(1/(moment_inertia_zz*e_modulus))*(r_3 * ( x^3/ 6 -x_3) + r_2 * ( x^3/ 6 -x_2) + r_1*(x^3/6 - x_1)- (q*x^4)/24)
        
    return v
x = 0
deflection_span = []
while x<= aileron_span:
    deflection_span.append(deflection(x))

print(deflection_span)
    
    


    
    