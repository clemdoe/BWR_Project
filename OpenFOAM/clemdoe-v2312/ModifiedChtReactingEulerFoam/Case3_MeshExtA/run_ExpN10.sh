#!/bin/bash
#SBATCH --job-name=ExpN11_3
##SBATCH --ntasks=5
#SBATCH --nodes=1           ## otra opcion es pedir una cantidad de nodos
#SBATCH --tasks-per-node=2    ## ver si conviene lanzar 1 o 4 procesos por nodo
#SBATCH --output=log
#SBATCH --error=errorLOG
#SBATCH --mail-user=dmgodino@gmail.com   --mail-type=end   # envia un mail al finalizar


#mpirun chtReactingTwoPhaseEulerFoam_Modif1 -parallel  ##buoyantPimpleFoam
mpirun chtReactingTwoPhaseEulerFoam_Modif1 -parallel ##buoyantPimpleFoam

