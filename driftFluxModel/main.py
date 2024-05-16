#Drift flux model for BWR fuel assembly class
# Used to calculate the drift velocity of the mixture in the fuel assembly
# Authors : Clément Huet

from class_MVF import FVM
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

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
def getAreas(A, Phi2Phi, f, D_h, K_loss, DV, Dz):
    A_chap = A +  (Phi2Phi/2) * ((f / D_h) + (K_loss / Dz)) * DV
    return A_chap

#Function to calculate the thermodynamic quality of the mixture
def getThermodynamicQuality(U, H, P):
    T = getTemperature(P, H)
    gas = IAPWS97(T = T, x = 0)
    liquid = IAPWS97(T = T, x = 1)
    h_g_sat = gas.h
    h_l_sat = liquid.h
    h_m = IAPWS97(P = P, h = H).h
    x_th = h_m - h_l_sat / h_g_sat - h_l_sat #h_l_sat, h_g_sat are the specific enthalpies of the liquid and gas phases at saturation temperature
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
def getD_h(L,l,geoType,cladRadius,Phi):
    if geoType == "square":
        return ((l*L)-(np.pi*cladRadius**2))/(cladRadius*2*np.pi)
    if geoType =="cylinder":
        print('Cylinder geometry not implemented yet')
        return

## Parameters of the system
# Equation resolution parameters
eps = 10**(-6)
N_iterations = 10000

# Constant of the problem
N_vol = 10 #Number of volumes for the discretization using the FVM class
Phi = 1 #Porosity
H = 2 #Height of the fuel rod m
l = 14.04*10**(-3) #Side of the square fuel rod m
L = 14.04*10**(-3) #Side of the square fuel rod m
Phi2Phi = 0.5 #Two-phase friction multiplier ???????,
f = 0.5 #correction factor for drag coefficient ???????
K_loss = 0.5 #loss coefficient ???????
g = 9.81 #gravity m/s2
q__ = 3000000 #volumetric heat generation rate W/m3
cladRadius = 6.52*10**(-3) #External radius of the clad m

#Initial/boundary conditions of the system
rho_l_start= 1000 #kg/m3
rho_g_start = 1 #kg/m3
epsilon_start = 0.00001
U_start = 7 #m/2
T_inlet = 500 #K
P_oulet = 10.800000 #MPa
P_inlet = 10.800000 #MPa
h_start = IAPWS97(T = T_inlet, P = P_inlet).h #J/kg

#Calulated values
DV = (H/N_vol)*((l*L)-(np.pi*cladRadius**2)) #Volume of the control volume m3
Area = ((l*L)-(np.pi*cladRadius**2)) #Area of the control volume m2
D_h = getD_h(L,l,"square",cladRadius,Phi) #Hydraulic diameter m2
Dz = H/N_vol #Height of the control volume m

V_gj_start = getDriftVelocity(rho_g_start, rho_l_start, g, D_h) #m/s
C0_start = getC0(rho_g_start, rho_l_start)
Vgj_prime_start = V_gj_start + (C0_start -1) * U_start #m/s
rho_start = rho_l_start * epsilon_start + rho_g_start * (1 - epsilon_start)
Dhfg_start = 1000 #J/kg specific enthalpy of vaporization

#Initial fields of the system
U = np.ones(N_vol)*U_start
P = np.ones(N_vol)
H = np.ones(N_vol)
rho_g_old = np.ones(N_vol)*rho_g_start
rho_l_old = np.ones(N_vol)*rho_l_start
rho_old = np.ones(N_vol)*rho_start
V_gj_old = np.ones(N_vol)*V_gj_start
epsilon_old = np.ones(N_vol)*epsilon_start
areaMatrix = np.ones(N_vol)*Area
areaMatrix_old_ = [getAreas(areaMatrix[i], Phi2Phi, f, D_h, K_loss, DV, Dz) for i in range(N_vol)]
Dhfg = np.ones(N_vol)*Dhfg_start

C0 = np.ones(N_vol)*C0_start
x_th = np.ones(N_vol)
T = np.ones(N_vol)


for i in range(N_iterations):

    U_old = U
    P_old = P
    H_old = H

    #Solving the system of equations of the mixture model
    #Solving the equation of mass concervation to find the velocity of the mixture
    mixtureVelocityClass = FVM(A00 = 1, A01 = 0, Am0 = rho_old[0], Am1 = rho_old[0], D0 = U_old[0], Dm1 = 0, N_vol = 10, H = 2)
    #Solving the equation of momentum concervation to find the pressure of the mixture
    mixturePressureClass = FVM(A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)
    #Solving the equation of energy concervation to find the enthalpy of the mixture
    mixtureEnthalpyClass = FVM(A00 = 1, A01 = 1, Am0 = 1, Am1 = 1, D0 = 1, Dm1 = 1, N_vol = 10, H = 1)

    #Filling inside of the resolution matrix for the velocity, pressure and enthalpy
    for i in range(mixtureVelocityClass.N_vol-1):
        mixtureVelocityClass.set_ADi(i, ci = - rho_old[i-1],
            ai = rho_old[i],
            bi = 0,
            di =  0)

        mixturePressureClass.set_ADi(i, ci = areaMatrix[i-1],
            ai = areaMatrix[i],
            bi = 0,
            di =  (((epsilon_old[i]/(1-epsilon_old[i]))* rho_l_old[i]*rho_g_old[i]*V_gj_old[i]^2*areaMatrix[i])/rho_old[i]) - (((epsilon_old[i-1]/(1-epsilon_old[i-1]))* rho_l_old[i-1]*rho_g_old[i-1]*V_gj_old[i-1]^2*areaMatrix[i-1])/rho_old[i-1]) - ((rho_old[i]- rho_old[i-1])* g * DV / 2) - (rho_g_old[i] * U_old[i] * areaMatrix_old_[i]) - (rho_g_old[i-1] * U_old[i-1] * areaMatrix_old_[i-1])              )
        
    #Solving the system of equations for the velocity and pressure
    mixtureVelocityClass.verticalResolution()
    U = mixtureVelocityClass.h
    mixturePressureClass.verticalResolution()
    P = mixturePressureClass.h

    #Filling inside of the resolution matrix for the enthalpy,
    #the mixture enthalpy needs to be solve after de U field because in contain U* and U
    for i in range(mixtureVelocityClass.N_vol-1):
        mixtureEnthalpyClass.set_ADi(i, ci = - rho_old[i-1] * U_old[i-1] * areaMatrix_old_[i-1] * U[i-1],
            ai = rho_old[i] * U_old[i] * areaMatrix_old_[i] * U[i],
            bi = 0,
            di =  q__ * DV - ((epsilon_old[i] * rho_l_old[i] * rho_g_old[i] * Dhfg[i] * V_gj_old[i])/rho_old[i]) - ((epsilon_old[i-1] * rho_l_old[i-1] * rho_g_old[i-1] * Dhfg[i-1] * V_gj_old[i-1])/rho_old[i-1]) + ((U_old[i] + epsilon_old[i]*(rho_l_old[i] - rho_g_old[i])*V_gj_old[i] / rho_old[i]) + (U_old[i-1] + epsilon_old[i-1]*(rho_l_old[i-1] - rho_g_old[i-1])*V_gj_old[i-1] / rho_old[i-1])) * (P_old[i] * areaMatrix[i] - P_old[i-1] * areaMatrix[i-1]) / 2 )


        #Update Areas, density, drift velocity and void fraction
        areaMatrix_old_ = getAreas(areaMatrix[i], Phi2Phi, f, D_h, K_loss, DV, Dz)
        C0[i] = getC0(rho_g_old[i], rho_l_old[i])
        x_th[i] = getThermodynamicQuality(U[i], H[i], P[i])
        V_gj_old[i] = getDriftVelocity(rho_g_old[i], rho_l_old[i], g, D_h)
        Vgj_prime[i] = Vgj + (C0 -1) * U_old[i]
        epsilon_old[i] = getVoidFraction(rho_g_old[i], rho_l_old[i], U[i], H[i], P[i], U_old[i], rho_old[i], D_h, g)
        T[i] = getTemperature(P[i], H[i])
        rho_g_old[i], rho_l_old[i], rho_old[i] = getDensity(U[i], H[i], P[i])

    mixtureEnthalpyClass.verticalResolution()
    H = mixtureEnthalpyClass.h

    if np.linalg.norm(U - U_old) < eps or np.linalg.norm(P - P_old) < eps or np.linalg.norm(H - H_old) < eps:
        print(f"Itération number: {i}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")
        print(U,P,H)
        print(mixtureVelocityClass.A)
        break

    elif i == N_iterations - 1:
        raise ValueError("The system did not converge")
    
    else:
        U_residual = np.linalg.norm(U - U_old)
        P_residual = np.linalg.norm(P - P_old)
        H_residual = np.linalg.norm(H - H_old)
        print(f"Itération number: {i}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")