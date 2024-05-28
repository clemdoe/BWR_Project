/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtReactingTwoPhaseEulerFoam

Description
    Solver for a system of 2 compressible fluid phases with a common pressure,
    but otherwise separate properties. The type of phase model is run time
    selectable and can optionally represent multiple species and in-phase
    reactions. The phase system is also run time selectable and can optionally
    represent different types of momentun, heat and mass transfer.

\*---------------------------------------------------------------------------*/

//Include propios del Dos-fluidos
#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "fvcSmooth.H"


// ----------------------------------------------------- //
//Include propios del CHT
#include "regionProperties.H"
#include "solidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "fixedGradientFvPatchFields.H"
#include "compressibleCourantNoFluid.H"
#include "solidRegionDiffNo.H"  //Incluido para el solido
#include "fvOptions.H"
#include "coordinateSystem.H"

// ----------------------------------------------------- //
// Controles de tiempo y pimple del dos fluidos y el CHT
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"

#include "pimpleMultiRegionControl.H"
#include "pressureControl.H"
// ----------------------------------------------------- //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "aclarDario.H"
    #include "createTime.H"

    // ---------------- Malla y campos ----------------------- //
    // Se cargan las malla y los campos de las diferentes regiones y del dos fluido
    Foam::regionProperties rp(runTime);

    // --- Malla y campo region fluido ---
    #include "createFluidMeshes.H"
    #include "createFluidFields.H"
    // --- Malla y campo region solido ---
    #include "createSolidMeshes.H"
    #include "createSolidFields.H"
    // --- Malla y campo region escalar (ESTE REGION OJO PORQUE NO LA USO)
    #include "createScalarMesh.H"
    #include "createScalarFields.H"
    // --- Malla y campo region dos fluidos ---
    #include "createMesh.H"
    //#include "createControl.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    pimpleControl pimple(mesh);
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    // ------------------------------------------------------ //


    // ------------------------------------------------------ //
    #include "initContinuityErrsFluid.H"
    #include "createFluidPressureControlsFluid.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    //#include "setInitialMultiRegionDeltaT.H"
    // ------------------------------------------------------ //


    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );

    //#include "pUf/createRDeltaTf.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        int nEnergyCorrectors
        (
            pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)
        );

        //La version original tenia incorporado el LTS
        //Que es algo de Local Time Step.
        //Ver si tiene sentido o no incorporarlo
        // ------------------------------------------------------ //
        // Se lee los diferentes parametros de tiempo y se elije el Courant maximo
        #include "readSolidTimeControls.H"
        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "CourantNos.H"
        #include "setDeltaTMultiRegion.H"
        // ------------------------------------------------------ //



        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            Info<< "\nSolving for two-fluid region " << endl;

            //Info<< "\nSolving for two-fluid region " << pimple.control_ << endl;

            fluid.solve();
            fluid.correct();

            #include "YEqns.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "EEqns.H"
                #include "pUf/pEqn.H"
            }
            else
            {
                #include "pU/UEqns.H"
                #include "EEqns.H"
                #include "pU/pEqn.H"
           }

            fluid.correctKinematics();

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }

            // ------------------------------------------------------ //
            Info << "\nResuelvo las regiones de acoplamiento " << endl;

            //Resuelvo las regiones del CHT - Tres regiones
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }

            //Esta region resuelvo la ecuacion de energia para una regions
            forAll(scalarRegions,i)
            {
                Info<< "\nSolving for scalar region "
                    << scalarRegions[i].name() << endl;
                #include "setRegionScalarFields.H"
                #include "scalarSolve.H"
            }
            // ------------------------------------------------------ //
        }

        // ------------------------------------------------------------ //
        // Actualizo los campo phiSteam y phiWater
        // Esto lo utilizo para calcular el caudal que sale por arriba
        phiSteam = fluid.phase1().phi()*linearInterpolate(fluid.phase1().thermoRef().rho()*fluid.phase1());

        phiWater = fluid.phase2().phi()*linearInterpolate(fluid.phase2().thermoRef().rho()*fluid.phase2());
        // ------------------------------------------------------------ //


        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
