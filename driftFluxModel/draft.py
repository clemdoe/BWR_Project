#Drift flux model for BWR fuel assembly class
# Used to calculate the drift velocity of the mixture in the fuel assembly
# Authors : Clément Huet

from class_MVF import FVM
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

rho_liquide_sat = IAPWS97(P=791.7e-3, x=0).rho
rho_vapeur_sat = IAPWS97(P=791.7e-3, x=1).rho
masse_tot = rho_liquide_sat*0.1 + rho_vapeur_sat*0.9

x = rho_vapeur_sat * 0.9 / masse_tot
print(x)
print(f'Masse totale : {masse_tot:.3f} kg.')