#Drift flux model for BWR fuel assembly class
# Used to calculate the drift velocity of the mixture in the fuel assembly
# Authors : Clement Huet

from class_MVF import MVF
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

## Parameters of the system
# Equation resolution parameters
eps = 10**-6
N_iterations = 10000

#Initial conditions of the system
rho_s = 2
rho_n = 2
epsilon = 0.0001
V_init = 12

# Meshing conditions
N_vol = 10

#Calculated parameters

U = np.ones(N_vol)
P = np.ones(N_vol)
H = np.ones(N_vol)


for i in range(N_iterations):

    U_old = U
    P_old = P
    H_old = H

    #Setting density, drift velocity, void fraction with consititutives equations
    rho = rho_s * (1 - U) + rho_n * U
    epsilon = 
    V_gj = 


    #Solving the system of equations of the mixture model
    #Solving the equation of mass concervation to find the velocity of the mixture
    mixtureVelocityClass = MVF(ai = rho_n, bi = 0, ci = rho_s, di = 0, A00 = 1, A01 = 0, Am0 = rho_s, Am1 = rho_s, D0 = V_init, Dm1 = 0, N_vol = 10, H = 2)
    mixtureVelocityClass.verticalResolution()
    U = mixtureVelocityClass.h
    print(U)
    
    #Solving the equation of momentum concervation to find the pressure of the mixture
    mixturePressureClass = MVF(ai = 1, bi = 1, ci = 1, di = 1, A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)
    mixturePressureClass.verticalResolution()
    P = mixturePressureClass.h

    #Solving the equation of energy concervation to find the enthalpy of the mixture
    mixtureEnthalpyClass = MVF(ai = 1, bi = 1, ci = 1, di = 1, A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)
    mixtureEnthalpyClass.verticalResolution()
    H = mixtureEnthalpyClass.h
    
    if np.linalg.norm(U - U_old) < eps or np.linalg.norm(P - P_old) < eps or np.linalg.norm(H - H_old) < eps:
        print(f"Itération number: {i}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")
        print(U,P,H)
        print(mixtureVelocityClass.A)
        break

    if i == N_iterations - 1:
        raise ValueError("The system did not converge")
    
    else:
        U_residual = np.linalg.norm(U - U_old)
        P_residual = np.linalg.norm(P - P_old)
        H_residual = np.linalg.norm(H - H_old)
        print(f"Itération number: {i}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")

    

