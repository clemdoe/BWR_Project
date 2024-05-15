#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os 
import numpy as np
import matplotlib.pyplot as plt

# Data Choice 
nomSimulation = "CHF_Temperature/"
fixedPower = False # True = fixedPower, False = fixedTemperature
tMin = 1
tMax = 77
tStep = 1

# Changing path to the simulation (to be CHANGED if needed)
#os.chdir("/home/gauthierlazare/Documents/GeN-Foam-develop/GauthierGeNFoam/"+nomSimulation)
print (os.getcwd())

# File names useful
cheminFichier = "/fluidRegion/"
alphaFile = "alpha.vapour"
saturationTemperatureFile = "T.interface"
heatFluxFile = "heatFlux.structure"
wallTemperatureFile = "T.fixedTemperature"

# Functions 
def imprimer(data):
	text = str(data)+"\n"
	return text

def ecrire(liste,nom):
	fichier = open(nom, "w")
	for data in liste:
		fichier.write(imprimer(data))
	fichier.close()


# Code 
N = (tMax-tMin)//tStep+1
t = [tMin+i*tStep for i in range(N)]

alpha = [0]*N
Tsat = [0]*N
Tw = [0]*N
hF = [0]*N

for k in range(N):
	ti = t[k]
	file = open(str(ti)+cheminFichier+alphaFile,"r")
	lines = file.readlines()
	M = len(lines)
	if (M>=35):
		alpha[k]=float(lines[35])
	else:
		alpha[k]=0
	file.close()
	file = open(str(ti)+cheminFichier+saturationTemperatureFile,"r")
	lines = file.readlines()
	M = len(lines)
	if (M>=35):
		Tsat[k]=float(lines[35])
	else:
		Tsat[k]=0
	file.close()

	file = open(str(ti)+cheminFichier+wallTemperatureFile,"r")
	lines = file.readlines()
	M = len(lines)
	if (M>=35):
		Tw[k]=float(lines[35])
	else:
		Tw[k]=0
	file.close()

	file = open(str(ti)+cheminFichier+heatFluxFile,"r")
	lines = file.readlines()
	M = len(lines)
	if (M>=35):
		hF[k]=float(lines[35])
	else:
		hF[k]=0
	file.close()

result = "Results/"

try:
    # Create target Directory
    os.mkdir(os.getcwd()+'/'+result)
    print("Directory " , result ,  " Created ") 
except FileExistsError:
    print("Directory " , result ,  " already exists")




ecrire(alpha,result+"alpha.txt")
ecrire(Tsat,result+"Tsat.txt")
ecrire(Tw,result+"Twall.txt")
ecrire(hF,result+"HeatFlux.txt")

Twnorm=np.array(Tw)-np.array(Tsat)
hfarray=np.array(hF)
plt.loglog(Twnorm, hfarray)
plt.show()