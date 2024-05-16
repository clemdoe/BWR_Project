#Drift flux model for BWR fuel assembly class
# Used to calculate the drift velocity of the mixture in the fuel assembly
# Authors : Clément Huet

from class_MVF import FVM
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

## Parameters of the system
# Equation resolution parameters
eps = 10**(-6)
N_iterations = 10000

#Initial conditions of the system
rho_l_s = 1000
rho_g_s = 1
epsilon = 0.0001
V_init = 12

# Constant of the problem
N_vol = 10 #Number of volumes for the discretization using the FVM class
Phi = 1 #Porosity
H = 2 #Height of the fuel rod
l = 0.5 #Side of the square fuel rod
L = 0.5 #Side of the square fuel rod
Phi2Phi = 0.5 #Two-phase friction multiplier
f = 0.5 #correction factor for drag coefficient
D_h = 0.5 #hydraulic-equivalent diameter
K_loss = 0.5 #loss coefficient
g = 9.81 #gravity
q__ = 0.5 #volumetric heat generation rate
cladRadius = 1 #External radius of the clad

#Calulated values
DV = (H/N_vol)*l*L #Volume of the control volume

#Initial fields of the system
U = np.ones(N_vol)
P = np.ones(N_vol)
H = np.ones(N_vol)
rho_g = np.ones(N_vol)
rho_l = np.ones(N_vol)
rho = np.ones(N_vol)
V_gj = np.ones(N_vol)
epsilon = np.ones(N_vol)

for i in range(N_iterations):

    U_old = U
    P_old = P
    H_old = H
    rho_g_old = rho_g
    rho_l_old = rho_l
    rho_old = rho
    V_gj_old = V_gj
    epsilon_old = epsilon

    rho_g[0] = 

    #Solving the system of equations of the mixture model
    #Solving the equation of mass concervation to find the velocity of the mixture
    mixtureVelocityClass = FVM(ai = rho_n, bi = 0, ci = rho_s, di = 0, A00 = 1, A01 = 0, Am0 = rho_s, Am1 = rho_s, D0 = V_init, Dm1 = 0, N_vol = 10, H = 2)
    #Solving the equation of momentum concervation to find the pressure of the mixture
    mixturePressureClass = FVM(ai = 1, bi = 1, ci = 1, di = 1, A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)
    #Solving the equation of energy concervation to find the enthalpy of the mixture
    mixtureEnthalpyClass = FVM(ai = 1, bi = 1, ci = 1, di = 1, A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)

    #Filling inside of the resolution matrix for the velocity, pressure and enthalpy
    for i in range(mixtureVelocityClass.N_vol-1):
        mixtureVelocityClass.set_ADi(i, ci = ,
            ai = ,
            bi = ,
            di =  )

        mixturePressureClass.set_ADi(i, ci = ,
            ai = ,
            bi = ,
            di =  )
        
        mixtureEnthalpyClass.set_ADi(i, ci = ,
            ai = ,
            bi = ,
            di =  )

        #Finding the density and temperature of the mixture
        T = getTemperature(P[i], H[i])
        rho_g[i], rho_l[i], rho[i] = getDensity(U[i], H[i], P[i])

        #Computing the areas for the velocity, pressure and enthalpy
        A_ = getAreas(A, Phi2Phi, f, D_h, K_loss, DV)
        Vgj[i] = getDriftVelocity(rho_g[i], rho_l[i], g, D_h)
        Vgj_prime[i] = Vgj + (C0 -1) * U_old[i]
        epsilonMatrix[i] = getVoidFraction(rho_g[i], rho_l[i], U[i], H[i], P[i], U_old[i], rho[i], D_h, g)

    mixtureVelocityClass.verticalResolution()
    U = mixtureVelocityClass.h
    mixturePressureClass.verticalResolution()
    P = mixturePressureClass.h
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

    

##Creation of function
#FUnction to calculate the temperature of the mixture
def getTemperature(P, H):
    T = IAPWS97(P = P, h = H).T
    return T

#Function to calculate the density of the mixture
def getDensity(U, H, P):
    T = getTemperature(P, H)
    gas = IAPWS97(T = T, x = 0)
    liquid = IAPWS97(T = T, x = 1)
    rho_g = gas.rho
    rho_l = liquid.rho
    rho = rho_s * (1 - U) + rho_n * U
    return rho_g, rho_l, rho

#Function to calculate the areas of the mixture
def getAreas(A, Phi2Phi, f, D_h, K_loss, DV):
    A_chap = A +  (Phi2Phi/2) * ((f / D_h) + (K_loss /Dz)) * DV
    return A_chap

#Function to calculate the thermodynamic quality of the mixture
def getThermodynamicQuality(U, H, P):
    T = getTemperature(P, H)
    gas = IAPWS97(T = T, x = 0)
    liquid = IAPWS97(T = T, x = 1)
    h_g_sat = gas.h
    h_l_sat = liquid.h
    h_m = IAPWS97(P = P, h = H).h
    x_th = h_m - h_l_sat / h_g_sat - h_l_sat
    return x_th

#Function to calculate the drift velocity of the mixture
def getDriftVelocity(rho_g, rho_l, g, D_h):
    return 0.188 * np.sqrt(((rho_l - rho_g) * g * D_h ) / rho_g )

#Function to calculate the constant C0 called the distribution parameter
def getC0(rho_g, rho_l):
    return 1.2 - 0.2*np.sqrt(rho_g / rho_l)

#Function to calculate the void fraction of the mixture
def getVoidFraction(rho_g, rho_l, U, H, P, U_old, rho, D_h, g):
    C0 = getC0(rho_g, rho_l)
    x_th = getThermodynamicQuality(U, H, P)
    Vgj = getDriftVelocity(rho_g, rho_l, g, D_h)
    Vgj_prime = Vgj + (C0 -1) * U_old
    epsilon = x_th / (C0 * (x_th + (rho_g/rho_l) * (1 - x_th)) + (rho_g * Vgj_prime / rho * U))
    return epsilon

#Function to calculate the hydraulic diameter
def getD_h():
    if geoType == "square":
        return (L*l*Phi)/(cladRadius*2*np.pi())
    if geoType =="cylinder":
        print('Cylinder geometry not implemented yet')
        break