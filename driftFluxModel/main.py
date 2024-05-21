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
def getDensity(epsilon, H, P):
    print(H,P)
    T = getTemperature(P, H)
    state = IAPWS97(T = T, x = epsilon)
    rho_g = state.Vapor.rho
    rho_l = state.Liquid.rho
    rho = rho_g * (1 - epsilon) + rho_l * epsilon
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
Height = 2 #Height of the fuel rod m
l = 14.04*10**(-3) #Side of the square fuel rod m
L = 14.04*10**(-3) #Side of the square fuel rod m
Phi2Phi = 0.5 #Two-phase friction multiplier ???????,
f = 0.5 #correction factor for drag coefficient ???????
K_loss = 0.5 #loss coefficient ???????
g = 9.81 #gravity m/s2
q__ = 300000 #volumetric heat generation rate W/m3
cladRadius = 6.52*10**(-3) #External radius of the clad m

#Initial/boundary conditions of the system
rho_l_start= 1000 #kg/m3
rho_g_start = 1 #kg/m3
epsilon_start = 0.00001
U_start = 7 #m/2
T_inlet = 500 #K
P_outlet = 10.800000 #MPa
P_inlet = 10.800000 #MPa
h_start = IAPWS97(T = T_inlet, P = P_inlet).h #J/kg

#Calulated values
DV = (Height/N_vol)*((l*L)-(np.pi*cladRadius**2)) #Volume of the control volume m3
Area = ((l*L)-(np.pi*cladRadius**2)) #Area of the control volume m2
D_h = getD_h(L,l,"square",cladRadius,Phi) #Hydraulic diameter m2
Dz = Height/N_vol #Height of the control volume m

#PATHS Values
inletFlowRate = 9.506 #m3/s
Power = 0
outletPressure = 7.18 #MPa
FlowArea = 0.010334 #m2
Height = 4.53*10**(-3) #m

#Calculated PATHS values
DV = (Height/N_vol) * FlowArea #Volume of the control volume m3
Area = FlowArea #Area of the control volume m2
D_h = Area / np.sqrt(Area)
Dz = Height/N_vol #Height of the control volume m
U_start = inletFlowRate * Area / rho_l_start #m/s

V_gj_start = getDriftVelocity(rho_g_start, rho_l_start, g, D_h) #m/s
C0_start = getC0(rho_g_start, rho_l_start)
#Vgj_prime_start = V_gj_start + (C0_start -1) * U_start #m/s
Vgj_prime_start = V_gj_start + (C0_start -1) * U_start #m/s
rho_start = rho_l_start * epsilon_start + rho_g_start * (1 - epsilon_start)
Dhfg_start = 500 #kJ/kg specific enthalpy of vaporization
h_inlet = IAPWS97(T = T_inlet, P = P_inlet).h #J/kg

#Initial fields of the system
U = np.ones(N_vol)#*U_start
P = np.ones(N_vol)
H = np.ones(N_vol)
rho_g_old = np.ones(N_vol)*rho_g_start
rho_l_old = np.ones(N_vol)*rho_l_start
rho_old = np.ones(N_vol)*rho_start
V_gj_old = np.ones(N_vol)*V_gj_start
Vgj_prime = np.ones(N_vol)*Vgj_prime_start
epsilon_old = np.ones(N_vol)*epsilon_start
areaMatrix = np.ones(N_vol)*Area
areaMatrix_old_ = [getAreas(areaMatrix[i], Phi2Phi, f, D_h, K_loss, DV, Dz) for i in range(N_vol)]
Dhfg = np.ones(N_vol)*Dhfg_start

C0 = np.ones(N_vol)*C0_start
x_th = np.ones(N_vol)
T = np.ones(N_vol)

for j in range(N_iterations):

    print(f"Itération number: {j}")

    U_old = U
    P_old = P
    H_old = H

    #Dm1 = ( q__ * DV - (epsilon_old[-1]*rho_l_old[-1]*rho_g_old[-1]*Dhfg[-1]*V_gj_old[-1]/rho_old[-1]) + (epsilon_old[-2]*rho_l_old[-2]*rho_g_old[-2]*Dhfg[-2]*V_gj_old[-2]/rho_old[-2]) + (1/2) * (P_old[-1]*areaMatrix[-1] - P_old[-2]*areaMatrix[-2]) * ( (U_old[-1] + epsilon_old[-1]*(rho_l_old[-1]-rho_g_old[-1])*V_gj_old[-1]/rho_old[-1]) + (U_old[-2] + epsilon_old[-2]*(rho_l_old[-2]-rho_g_old[-2])*V_gj_old[-2]/rho_old[-2])) )

    #Solving the system of equations of the mixture model
    #Solving the equation of mass concervation to find the velocity of the mixture
    mixtureVelocityClass = FVM(A00 = 1, A01 = 0, Am0 = - rho_old[0], Am1 = rho_old[1], D0 = U_start, Dm1 = 0, N_vol = 10, H = Height)
    #Solving the equation of momentum concervation to find the pressure of the mixture
    mixturePressureClass = FVM(A00 = - areaMatrix[0], A01 = areaMatrix[1], Am0 = 0, Am1 = 1, D0 = (((epsilon_old[1]/(1-epsilon_old[1]))* rho_l_old[1]*rho_g_old[1]*(V_gj_old[1]**2)*areaMatrix[1])/rho_old[1]) - (((epsilon_old[0]/(1-epsilon_old[0]))* rho_l_old[0]*rho_g_old[0]*(V_gj_old[0]**2)*areaMatrix[0])/rho_old[0]) - ((rho_old[1]- rho_old[0])* g * DV / 2) - (rho_g_old[1] * U_old[1] * areaMatrix_old_[1]) - (rho_g_old[0] * U_old[0] * areaMatrix_old_[0]), Dm1 = P_outlet, N_vol = 10, H = Height)
    
    #Filling inside of the resolution matrix for the velocity, pressure and enthalpy
    for i in range(1, mixtureVelocityClass.N_vol-1):
        mixtureVelocityClass.set_ADi(i, ci = - rho_old[i-1],
            ai = rho_old[i],
            bi = 0,
            di =  0)
        
    #Solving the system of equations for the velocity
    mixtureVelocityClass.verticalResolution()
    #print(f'mixtureVelocityClass.A: {mixtureVelocityClass.A}, mixtureVelocityClass.D: {mixtureVelocityClass.D}')
    U = mixtureVelocityClass.h

    for i in range(1, mixturePressureClass.N_vol-1):
        mixturePressureClass.set_ADi(i, ci = 0,
            ai = - areaMatrix[i-1],
            bi = areaMatrix[i],
            di =  (((epsilon_old[i]/(1-epsilon_old[i]))* rho_l_old[i]*rho_g_old[i]*(V_gj_old[i]**2)*areaMatrix[i])/rho_old[i]) - (((epsilon_old[i-1]/(1-epsilon_old[i-1]))* rho_l_old[i-1]*rho_g_old[i-1]*(V_gj_old[i-1]**2)*areaMatrix[i-1])/rho_old[i-1]) - ((rho_old[i]- rho_old[i-1])* g * DV / 2) - (rho_g_old[i] * U_old[i] * areaMatrix_old_[i] * U[i]) - (rho_g_old[i-1] * U_old[i-1] * areaMatrix_old_[i-1] * U[i-1] ))
        
    #Solving the system of equations for the pressure
    mixturePressureClass.verticalResolution()
    #print(f'mixturePressureClass.A: {mixturePressureClass.A}, mixturePressureClass.D: {mixturePressureClass.D}')
    P = mixturePressureClass.h

    #Solving the equation of energy concervation to find the enthalpy of the mixture
    i = -1
    Dm1 = q__ * DV - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]/rho_old[i-1]) + (1/2) * (P[i]*areaMatrix[i] - P[i-1]*areaMatrix[i-1]) * ( (U[i] + epsilon_old[i]*(rho_l_old[i]-rho_g_old[i])*V_gj_old[i]/rho_old[i]) + (U[i-1] + epsilon_old[i-1]*(rho_l_old[i-1]-rho_g_old[i-1])*V_gj_old[i-1]/rho_old[i-1]) )
    #print(f'Dm1: {Dm1}, epsilon_old[i]: {epsilon_old[i]}, rho_l_old[i]: {rho_l_old[i]}, rho_g_old[i]: {rho_g_old[i]}, Dhfg[i]: {Dhfg[i]}, V_gj_old[i]: {V_gj_old[i]}, rho_old[i]: {rho_old[i]}, P[i]: {P[i]}, areaMatrix[i]: {areaMatrix[i]}, U[i]: {U[i]}, U[i-1]: {U[i-1]}')
    mixtureEnthalpyClass = FVM(A00 = 1, A01 = 0, Am0 = - rho_old[-2] * U[-2] * areaMatrix_old_[-2], Am1 = rho_old[-1] * U[-1] * areaMatrix_old_[-1], D0 = h_inlet, Dm1 = Dm1, N_vol = 10, H = 1)
    #Filling inside of the resolution matrix for the enthalpy,
    #the mixture enthalpy needs to be solve after de U field because in contain U* and U
    print(f'epsilon_old: {epsilon_old},\n rho_g_old: {rho_g_old}, \n rho_l_old: {rho_l_old},\n  rho_old: {rho_old},\n V_gj_old: {V_gj_old}, \n Vgj_prime: {Vgj_prime},\n areaMatrix: {areaMatrix},\n Dhfg: {Dhfg},\n C0: {C0}, \n x_th: {x_th},\n T: {T}, P: {P}, U: {U}, H: {H}    ')

    for i in range(1, mixtureEnthalpyClass.N_vol-1):

        #DI = q__ * DV - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]/rho_old[i-1]) + (1/2) * (P[i]*areaMatrix[i] - P[i-1]*areaMatrix[i-1]) * ( (U[i] + epsilon_old[i]*(rho_l_old[i]-rho_g_old[i])*V_gj_old[i]/rho_old[i]) + (U[i-1] + epsilon_old[i-1]*(rho_l_old[i-1]-rho_g_old[i-1])*V_gj_old[i-1]/rho_old[i-1]) )
        DI = (P[i]*areaMatrix[i] - P[i-1]*areaMatrix[i-1]) #* ( (U[i] + epsilon_old[i]*(rho_l_old[i]-rho_g_old[i])*V_gj_old[i]/rho_old[i]) + (U[i-1] + epsilon_old[i-1]*(rho_l_old[i-1]-rho_g_old[i-1])*V_gj_old[i-1]/rho_old[i-1]) )
    
        print(f'di: {DI},  P[i]: {P[i]},  P[i-1]: {P[i-1]}, areaMatrix[i]: {areaMatrix[i]},  areaMatrix[i-1]: {areaMatrix[i-1]}, U[i]: {U[i]}, U[i-1]: {U[i-1]}')
        #print( f'di: {di}, epsilon_old[i]: {epsilon_old[i]}, rho_l_old[i]: {rho_l_old[i]}, rho_g_old[i]: {rho_g_old[i]}, Dhfg[i]: {Dhfg[i]}, V_gj_old[i]: {V_gj_old[i]}, rho_old[i]: {rho_old[i]}, epsilon_old[i-1]: {epsilon_old[i-1]}, rho_l_old[i-1]: {rho_l_old[i-1]}, rho_g_old[i-1]: {rho_g_old[i-1]}, Dhfg[i-1]: {Dhfg[i-1]}, V_gj_old[i-1]: {V_gj_old[i-1]}, rho_old[i-1]: {rho_old[i-1]}, P[i]: {P[i]}, P[i-1]: {P[i-1]}, areaMatrix[i]: {areaMatrix[i]},  areaMatrix[i-1]: {areaMatrix[i-1]}, U[i]: {U[i]}, U[i-1]: {U[i-1]}')
        mixtureEnthalpyClass.set_ADi(i, ci = - rho_old[i-1] * U[i-1] * areaMatrix_old_[i-1],
            ai = rho_old[i] * U[i] * areaMatrix_old_[i],
            bi = 0,
            di =  DI)
        
    mixtureEnthalpyClass.verticalResolution()
    print(f'mixtureEnthalpyClass.A: {mixtureEnthalpyClass.A}, mixtureEnthalpyClass.D: {mixtureEnthalpyClass.D}')
    #H = [mixtureEnthalpyClass.h[i]*rho_old[i]*0.001 for i in range(N_vol)]
    H = [mixtureEnthalpyClass.h[i] for i in range(N_vol)]
    #print(f'U : {mixtureVelocityClass.A}, P: {mixturePressureClass.A}, H: {mixtureEnthalpyClass.A}')
    print(f'U: {U}, P: {P}, H: {H}')
    for i in range(mixturePressureClass.N_vol):
        print(f'Mise à jour des variables constitutives. i = {i}')
        rho_g_old[i], rho_l_old[i], rho_old[i] = getDensity(epsilon_old[i], H[i], P[i])
        V_gj_old[i] = getDriftVelocity(rho_g_old[i], rho_l_old[i], g, D_h)
        C0[i] = getC0(rho_g_old[i], rho_l_old[i])
        Vgj_prime[i] = V_gj_old[i] + (C0[i] -1) * U[i]
        epsilon_old[i] = getVoidFraction(rho_g_old[i], rho_l_old[i], U[i], H[i], P[i], U_old[i], rho_old[i], D_h, g)
        areaMatrix_old_[i] = getAreas(areaMatrix[i], Phi2Phi, f, D_h, K_loss, DV, Dz)
        Dhfg[i] = 1000 #J/kg specific enthalpy of vaporization

    U_residual = np.linalg.norm(U - U_old)
    P_residual = np.linalg.norm(P - P_old)
    H_residual = np.linalg.norm(H - H_old)

    if np.linalg.norm(U - U_old) < eps or np.linalg.norm(P - P_old) < eps or np.linalg.norm(H - H_old) < eps:
        print(f"Itération number: {j}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")
        break

    elif j == N_iterations - 1:
        raise ValueError("The system did not converge")
    
    else:
        print(f"Itération number: {j}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")