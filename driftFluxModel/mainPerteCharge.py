#Drift flux model for BWR fuel assembly class
# Used to calculate the drift velocity of the mixture in the fuel assembly
# Authors : Clément Huet

from class_MVF import FVM
import numpy as np
from iapws import IAPWS97
import csv
import matplotlib.pyplot as plt

##Creation of function
#FUnction to calculate the temperature of the mixture
def getTemperature(P, H):
    T = IAPWS97(P = P, h = H).T
    return T

""" #Function to calculate the density of the mixture
def getDensity(epsilon, H, P):
    T = getTemperature(P, H)
    #print(f'Inside getDensity: epsilon: {epsilon}, H: {H}, P: {P}, T: {T}')
    rho_g_0 = 916.8 #kg/m3
    rho_l_0 = 1000 #kg/m3
    rho_g = rho_g_0 - 0.139 * (T - 273.15)
    rho_l = rho_l_0 - 0.4583 * (T - 273.15)
    rho = rho_l * (1 - epsilon) + rho_g * epsilon
    return rho_g, rho_l, rho """

#Function to calculate the density of the mixture
def getDensity(epsilon, H, P):
    T = getTemperature(P, H)
    #print(f'Inside getDensity: epsilon: {epsilon}, H: {H}, P: {P}, T: {T}')
    state = IAPWS97(T = T, x = epsilon)
    rho_g = state.Vapor.rho
    rho_l = state.Liquid.rho
    if rho_g == None:
        rho_g = 0
    rho = rho_l * (1 - epsilon) + rho_g * epsilon
    return rho_g, rho_l, rho

#Function to calculate the areas of the mixture
def getAreas(A, Phi2Phi, f, D_h, K_loss, DV, Dz):
    A_chap = A +  (Phi2Phi/2) * ((f / D_h) + (K_loss / Dz)) * DV
    return A_chap

#Function to calculate the thermodynamic quality of the mixture
def getThermodynamicQuality(H, P):
    T = getTemperature(P, H)
    """
    temp_sat_gas = IAPWS97(P = P, x = 0).T
    temp_sat_liquid = IAPWS97(P = P, x = 1).T
    h_g_sat = IAPWS97(P = P, T = temp_sat_gas).h
    h_l_sat = IAPWS97(P = P, T = temp_sat_liquid).h"""
    gas = IAPWS97(T = T, x = 0) 
    liquid = IAPWS97(T = T, x = 1)
    h_g_sat = gas.h
    h_l_sat = liquid.h
    print(f'H: {H}, P: {P}, T: {T}, h_l_sat: {h_l_sat}, h_g_sat: {h_g_sat}')
    h_m = IAPWS97(P = P, h = H).h
    print(f'h_m: {h_m}')
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
    x_th = getThermodynamicQuality(H, P)
    Vgj = getDriftVelocity(rho_g, rho_l, g, D_h)
    Vgj_prime = Vgj + (C0 -1) * U_old
    epsilon = x_th / (C0 * (x_th + (rho_g/rho_l) * (1 - x_th)) + (rho_g * Vgj_prime / rho * U))
    #print(f'Inside getVoidFraction: rho_g: {rho_g}, rho_l: {rho_l}, U: {U}, H: {H}, P: {P}, U_old: {U_old}, rho: {rho}, D_h: {D_h}, g: {g}, C0: {C0}, x_th: {x_th}, Vgj: {Vgj}, Vgj_prime: {Vgj_prime}, epsilon: {epsilon}')
    return epsilon

def get_parameters(P,H,rho_l,rho_g,rho,epsilon, D_h, g, U, U_old):
    rho_l_old = rho_l
    rho_g_old = rho_g
    rho_old = rho
    epsilon_old = epsilon
    #print(f'Inside get parameters: rho_l_old: {rho_l_old}, rho_g_old: {rho_g_old}, rho_old: {rho_old}, epsilon_old: {epsilon_old}')
    for i in range(1000):
        epsilon = getVoidFraction(rho_g_old, rho_l_old, U, H, P, U_old, rho_old, D_h, g)
        #print(f'Inside get parameters inside the boucle: rho_l_old: {rho_l_old}, rho_g_old: {rho_g_old}, rho_old: {rho_old}, epsilon: {epsilon}')  
        rho_g, rho_l, rho = getDensity(epsilon, H, P)
        if abs(rho_g - rho_g_old) < 10**(-3) and abs(rho_l - rho_l_old) < 10**(-3) and abs(rho - rho_old) < 10**(-3) and abs(epsilon - epsilon_old) < 10**(-3):
            break
        else:
            rho_g_old = rho_g
            rho_l_old = rho_l
            rho_old = rho
            epsilon_old = epsilon
        
    x_th = getThermodynamicQuality(H, P)
    C0 = getC0(rho_g, rho_l)
    Vgj = getDriftVelocity(rho_g, rho_l, g, D_h)
    
    return rho_g, rho_l, rho, epsilon, x_th, C0, Vgj

#Function to calculate the hydraulic diameter
def getD_h(L,l,geoType,cladRadius):
    if geoType == "square":
        return ((l*L)-(np.pi*cladRadius**2))/(cladRadius*2*np.pi)
    if geoType =="cylinder":
        print('Cylinder geometry not implemented yet')
        return
    
#Function to merge the 3 component of the matrix U,P,H
def createVar(U,P,H):
    var = np.zeros(3*len(U))
    for i in range(len(var)):
        if i < len(U):
            var[i] = U[i]
        elif i < 2*len(U):
            var[i] = P[i-len(U)]
        else:
            var[i] = H[i-2*len(U)]
    return var

#Function to split the 3 component of the matrix U,P,H
def splitVar(var):
    U = var[:len(var)//3]
    P = var[len(var)//3:2*len(var)//3]
    H = var[2*len(var)//3:]
    return list(U),list(P),list(H)

## Parameters of the system
# Equation resolution parameters
eps = 10**(-3)
N_iterations = 1000

# Constant of the problem
sizeMesh = 4 #Number of volumes for the discretization using the FVM class
N_vol = 12 #Number of volumes for the discretization using the FVM class
Phi = 1 #Porosity
Height = 2 #Height of the fuel rod m
l = 14.04*10**(-3) #Side of the square fuel rod m
L = 14.04*10**(-3) #Side of the square fuel rod m
Phi2Phi = 0.5 #Two-phase friction multiplier ???????
Phi_start = 50 #Lockart Martinelli parameter
f = 0.005 #correction factor for drag coefficient Lockart Martinelli
K_loss = 0 #loss coefficient 
g = 9.81 #gravity m/s2
q__ = 0 #volumetric heat generation rate W/m3
cladRadius = 6.52*10**(-3) #External radius of the clad m
waterGap = 0.5*10**(-3) #Gap between the clad and the water m
waterRadius =  cladRadius + waterGap #External radius of the water m
massFlowRate = 7000 #kg/m2/s

#Initial/boundary conditions of the system
rho_l_start= 1000 #kg/m3
rho_g_start = 916.8 #kg/m3
epsilon_start = 0
U_start = 7 #m/2
T_inlet = 500 #K
P_outlet = 10.800000*10**(6) #Pa
P_inlet = 10.800000*10**(6) #Pa
h_start = IAPWS97(T = T_inlet, P = P_inlet*10**(-6)).h *1000 #J/kg

#Calulated values
Area = ((np.pi*waterRadius**2)-(np.pi*cladRadius**2)) #Area of the control volume m2
DV = (Height/N_vol)*Area #Volume of the control volume m3
U_start = massFlowRate / (rho_l_start) #m/s
#Area = 2
D_h = Area / (np.pi*cladRadius) #Hydraulic diameter m2
#D_h = getD_h(L,l,"square",cladRadius,Phi) #Hydraulic diameter m2
Dz = Height/N_vol #Height of the control volume m

#print(f"DV: {DV}, Area: {Area}, D_h: {D_h}, Dz: {Dz}, Height: {Height}, L: {L}, l: {l}, rho_l_start: {rho_l_start}, rho_g_start: {rho_g_start}, epsilon_start: {epsilon_start}, U_start: {U_start}, T_inlet: {T_inlet}, P_outlet: {P_outlet}, P_inlet: {P_inlet}, h_start: {h_start}")

""" #PATHS Values
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
U_start = inletFlowRate * Area / rho_l_start #m/s """

V_gj_start = getDriftVelocity(rho_g_start, rho_l_start, g, D_h) #m/s
C0_start = getC0(rho_g_start, rho_l_start)
#Vgj_prime_start = V_gj_start + (C0_start -1) * U_start #m/s
Vgj_prime_start = V_gj_start + (C0_start -1) * U_start #m/s
rho_start = rho_g_start * epsilon_start + rho_l_start * (1 - epsilon_start)
Dhfg_start = 500 #kJ/kg specific enthalpy of vaporization
h_inlet = IAPWS97(T = T_inlet, P = P_inlet*10**(-6)).h * 1000 #J/kg

#Initial fields of the system
U = np.ones(sizeMesh)
P = np.ones(sizeMesh)
H = np.ones(sizeMesh)
rho_old = np.ones(sizeMesh)*rho_start
rho_g_old = np.ones(sizeMesh)*rho_g_start
rho_l_old = np.ones(sizeMesh)*rho_l_start
V_gj_old = np.ones(sizeMesh)*V_gj_start
Vgj_prime = np.ones(sizeMesh)*Vgj_prime_start
epsilon_old = np.ones(sizeMesh)*epsilon_start
areaMatrix = np.ones(sizeMesh)*Area
Dhfg = np.ones(sizeMesh)*Dhfg_start
Phi = np.ones(sizeMesh)*Phi_start
Friction = np.ones(sizeMesh)*f
areaMatrix_old_ = [getAreas(areaMatrix[i], Phi[i], Friction[i], D_h, K_loss, DV, Dz) for i in range(sizeMesh)]

C0 = np.ones(sizeMesh)*C0_start
x_th = np.ones(sizeMesh)
T = np.ones(sizeMesh)

for j in range(N_iterations):

    print(f"\n Begin itération number: {j}")
    
    U_old = U
    P_old = P
    H_old = H

    VAR_old = createVar(U_old,P_old,H_old)
    rho_old = createVar(rho_old, rho_old, rho_old)
    areaMatrix_old_ = createVar(areaMatrix_old_, areaMatrix_old_, areaMatrix_old_)
    areaMatrix = createVar(areaMatrix, areaMatrix, areaMatrix)
    print(f'AreaMatrix: {areaMatrix}', f'AreaMatrix_old_: {areaMatrix_old_}')

    DM1 = (1/2) * (P[-1]*areaMatrix[-1] - P[-2]*areaMatrix[-2]) * ( (U_old[-1]) + (U_old[-2]) )
    VAR_VFM_Class = FVM(A00 = 1, A01 = 0, Am0 = - rho_old[-2] * VAR_old[sizeMesh-2] * areaMatrix[-2], Am1 = rho_old[-1] * VAR_old[sizeMesh-1] * areaMatrix[-1], D0 = U_start, Dm1 = DM1, N_vol = 12, H = Height)
    VAR_VFM_Class.boundaryFilling()
    for i in range(1, VAR_VFM_Class.N_vol-1):
        #Inside the velocity submatrix
        if i < sizeMesh-1:
            VAR_VFM_Class.set_ADi(i, ci = - rho_old[i-1]*areaMatrix[i-1],
            ai = rho_old[i]*areaMatrix[i],
            bi = 0,
            di =  0)
        elif i == sizeMesh-1:
            VAR_VFM_Class.set_ADi(i, 
            ci = - rho_old[i-1]*areaMatrix[i-1],
            ai = rho_old[i]*areaMatrix[i],
            bi = 0,
            di =  0)

        #Inside the pressure submatrix
        elif i == sizeMesh:
            #DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 + areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 + areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
            VAR_VFM_Class.set_ADi(sizeMesh, 
            ci = 0,
            ai = - areaMatrix[i],
            bi = areaMatrix[i+1],
            di = - ((rho_old[i+1]- rho_old[i])* g * DV / 2))
        
            VAR_VFM_Class.fillingOutsideBoundary(i, i-sizeMesh,
            ai = - rho_old[i]*VAR_old[i-sizeMesh]*areaMatrix_old_[i],
            bi = rho_old[i+1]*VAR_old[i-sizeMesh]*areaMatrix_old_[i+1])

        elif i > sizeMesh and i < 2*sizeMesh-1:
            #DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 + areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 + areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
            VAR_VFM_Class.set_ADi(i, ci = 0,
            ai = - areaMatrix[i],
            bi = areaMatrix[i+1],
            di = - ((rho_old[i+1]- rho_old[i])* g * DV / 2))
        
            VAR_VFM_Class.fillingOutsideBoundary(i, i-sizeMesh,
            ai = - rho_old[i]*VAR_old[i-sizeMesh]*areaMatrix_old_[i],
            bi = rho_old[i+1]*VAR_old[i+1-sizeMesh]*areaMatrix_old_[i+1])

        elif i == 2*sizeMesh -1:
            VAR_VFM_Class.set_ADi(i, 
            ci = 0,
            ai = 1,
            bi = 0,
            di =  P_outlet)

            VAR_VFM_Class.fillingOutsideBoundary(2*sizeMesh -1, 2*sizeMesh -1 - sizeMesh,
            ai = 0,
            bi = 0)

        #Inside the enthalpy submatrix
        elif i == 2*sizeMesh:
            VAR_VFM_Class.set_ADi(2*sizeMesh, 
            ci = 0,
            ai = 1,
            bi = 0,
            di =  h_inlet)

        elif i > 2*sizeMesh and i < 3*sizeMesh:
            DI = (1/2) * (VAR_old[i-sizeMesh]*areaMatrix[i] - VAR_old[i-1-sizeMesh]*areaMatrix[i-1]) * ( (VAR_old[i-2*sizeMesh] ) + (VAR_old[i-1-2*sizeMesh]) )
            VAR_VFM_Class.set_ADi(i, ci =  - rho_old[i-1] * VAR_old[i-1-sizeMesh] * areaMatrix[i-1],
            ai = rho_old[i] * VAR_old[sizeMesh] * areaMatrix[i],
            bi = 0,
            di =  DI)

    print(f'VAR_VFM_Class.A: {VAR_VFM_Class.A}, VAR_VFM_Class.D: {VAR_VFM_Class.D}')
    VAR = VAR_VFM_Class.resoudre_h()
    U, P, H = splitVar(VAR)
    print(f'U: {U}, P: {P}, H: {H}')

    rho_old, rho_old0, rho_old0 = splitVar(rho_old)
    areaMatrix_old_, areaMatrix_old0, areaMatrix_old0 = splitVar(areaMatrix_old_)
    areaMatrix, areaMatrix0, areaMatrix0 = splitVar(areaMatrix)

    for i in range(sizeMesh):
        areaMatrix_old_[i] = getAreas(areaMatrix[i], Phi2Phi, f, D_h, K_loss, DV, Dz)
    
    for i in range(sizeMesh):
        rho_g_old[i], rho_l_old[i], rho_old[i], epsilon_old[i], x_th[i], C0[i], V_gj_old[i] = get_parameters(P[i]*10**(-6), H[i]*10**(-3), rho_l_old[i], rho_g_old[i], rho_old[i], epsilon_old[i], D_h, g, U[i], U_old[i])
        print(f"rho_g_old[i]: {rho_g_old[i]}, rho_l_old[i]: {rho_l_old[i]}, rho_old[i]: {rho_old[i]}, epsilon_old[i]: {epsilon_old[i]}, x_th[i]: {x_th[i]}, C0[i]: {C0[i]}, V_gj_old[i]: {V_gj_old[i]}")
    
    
    U_residual = np.linalg.norm(np.array(U) - np.array(U_old))
    P_residual = np.linalg.norm(np.array(P) - np.array(P_old))
    H_residual = np.linalg.norm(np.array(H) - np.array(H_old))

    if (np.linalg.norm(np.array(U) - np.array(U_old)) < eps) and (np.linalg.norm(np.array(P) - np.array(P_old)) < eps) and (np.linalg.norm(np.array(H) - np.array(H_old)) < eps):
        print(f"Itération number: {j}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")
        print(f"Convergence reached at iteration {j}")
        print(f'U: {U}, P: {P}, H: {H}')
        break

    elif j == N_iterations - 1:
        print(f'rho_old: {rho_old}, Area : {Area})')
        raise ValueError("The system did not converge")
        
    else:
        print(f"Itération number: {j}, U_residual: {U_residual}, P_residual: {P_residual}, H_residual: {H_residual}")
        print(f"Convergence not reached yet at iteration {j}")
