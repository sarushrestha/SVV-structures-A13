# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 10:37:54 2019

@author: Group A13
"""

import numpy as np

def reaction_forces(Ca, la, x1, x2, x3, xa, h, d1, d3, theta, P, q, E, I):
    """Calculates all reaction forces 
    based on geometery inputs and external forces"""
    
    equation_matrix = np.array([[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                               [1, 0, 0, 1, 0, 1, 0, np.sin(theta), 0, 0, 0, 0, (P*np.sin(theta)+q*la*np.cos(theta))], 
                               [0, 1, 0, 0, 1, 0, 1, np.cos(theta), 0, 0, 0, 0, (P*np.cos(theta)-q*la*np.sin(theta))],
                               
                               [-(Ca/4-h/2), 0, 0, -(Ca/4-h/2) ,0 , -(Ca/4-h/2), 0, (np.cos(theta)*h/2-np.sin(theta)*Ca/4), 0, 0, 0, 0, (P*np.cos(theta)*h/2*-P*np.sin(theta)*Ca/4)], 
                               [0, (x2-x1), 0, 0, 0, 0, -(x3-x2), (np.cos(theta)*xa/2), 0, 0, 0, 0, (-P*np.cos(theta)*xa/2+q*la*np.sin(theta)*(la/2-x2))], 
                               [-(x2-x1), 0, 0, 0, 0, (x3-x2), 0, -np.sin(theta)*xa/2, 0, 0, 0, 0, (P*np.sin(theta)*xa/2+q*la*np.cos(theta)*(la/2-x2))], 
                               
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x1, 1, -q*np.sin(theta)*((x1**4)/24)],  
                               [0, ((x2-x1)**3)/6, 0, 0, 0, 0, 0, ((np.cos(theta))*((xa/2)**3)/6), 0, 0, x2, 1, (-q*np.sin(theta)*((x2**4)/24))], 
                               [0, ((x3-x1)**3)/6, 0, 0, ((x3-x2)**3)/6, 0, 0, ((np.cos(theta))*((x3-x2+xa/2)**3)/6), 0, 0, x3, 1, (-q*np.sin(theta)*((x3**4)/24)+P*(np.cos(theta))*(x3-x2-xa/2)**3/6)], 
                               [0, 0, 0, 0, 0, 0, 0, 0, x1, 1, 0, 0, (-E*I*d1*+q*np.cos(theta)*(x1**4)/24)], 
                               [(((x2-x1)**3)/6), 0, 0, 0, 0, 0, 0, ((-np.sin(theta))*((xa/2)**3)/6), x2, 1, 0, 0, (q*np.cos(theta)*(x2**4)/24)], 
                               [(((x3-x1)**3)/6),0,0,(((x3-x2)**3)/6),0,0,0,((-np.sin(theta))*((x3-x2+xa/2)**3)/6),x3,1,0,0,(-E*I*d3*+q*np.cos(theta)*((x3**4)/24)+P/6*np.sin(theta)*(x3-x2-xa/2)**3)]])
    
    
    unknown_matrix = equation_matrix[:,:-1]
    constant_matrix = equation_matrix[:,-1]
    
    
    solution_matrix = np.linalg.solve(unknown_matrix,constant_matrix)
    
    solution_matrix = solution_matrix/1000
    
    (R1y, R1z, R2x, R2y, R2z, R3y, R3z, RI, c1, c2, c3, c4) = tuple(solution_matrix)
    
    print((R1y, R1z, R2x, R2y, R2z, R3y, R3z, RI, c1, c2, c3, c4))



(Ca, la, x1, x2, x3, xa, h, d1, d3, theta, P, q, E, t) = (0.515, 2.691, 0.174, 1.051, 2.512, 0.30, 0.248, 0.01034, 0.02066, 25*3.14159/180, 20600, 1000, 71*10**9, 0.0011)
I = (Ca*h**3)/12-((Ca-2*t)*(h-2*t)**3)/12
print(I)
reaction_forces(Ca, la, x1, x2, x3, xa, h, d1, d3, theta, P, q, E, I)