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
    scalarSolver.H resuelve dentro del reactingTwoPhaseEulerFoam el transporte de un escalar pasivo

Description
    Resuelvo dentro del reactingTwoPhaseEulerFoam la ecuacion de energia.
    OJO PORQUE NO ESTA BIEN IMPLEMENTADO ESTO!!!
    SI ESTA PREPARADO PARA PODER IMPLEMENTAR UN DETERMINADA ECUACION PERO
    LA QUE ESTA IMPEMENTADA NO FUNCIONA

\*---------------------------------------------------------------------------*/

//Info << "Datos del caso scalar" << endl;
//Info << "La energia cinetica: " << Kscal << endl;
//Info << "La densidad: " << rhoscal << endl;
//Info << "La gravedad: " << gscal << endl;
//Info << "El valor de h: " << hscal << endl;
//Info << "El valor de la velocidad: " << Uscal << endl;
//Info << "El valor de phi: " << phiscal << endl;
//Info << "El valor de Kappa: " << thermoScal.kappa() << endl;

 //fvm::ddt(rho, he) + fvm::div(phi, he)

Info << "El campo hscal inicial - max: " << max(hscal.primitiveField()) << "- min: " << min(hscal.primitiveField()) << endl;

    //while (pimple.correctNonOrthogonal())
   // {
          fvScalarMatrix hTEqn
          (
                      fvm::ddt(rhoscal, hscal) + fvm::div(phiscal, hscal,"div(phiscal,hscal)")
                    + fvc::ddt(rhoscal, Kscal) + fvc::div(phiscal, Kscal,"div(phiscal,Kscal)")
//                    + (
//                       hscal.name() == "e"
//                       ? fvc::div
//                       (
//                            fvc::absolute(phi/fvc::interpolate(rho), U),
//                            p,
//                            "div(phiv,p)"
//                       )
//                            : -dpdt
//                       )
                       - fvm::laplacian(thermoScal.kappa()/thermoScal.Cp(), hscal, "laplacian(kapparhoCp,hscal)")
                    ==
                       rhoscal*(Uscal&gscal)
                     + fvOptions(rhoscal, hscal)
                     );

                     //hTEqn.relax();

                     fvOptions.constrain(hTEqn);

                     hTEqn.solve();

                     fvOptions.correct(hscal);

//                     Info<< "Min/max T:" << min(thermoScal.T()).value() << ' '
//                         << max(thermoScal.T()).value() << endl;


                     thermoScal.correct();

//                     Info<< "Min/max T:" << min(thermoScal.T()).value() << ' '
//                         << max(thermoScal.T()).value() << endl;
   // }
                     Info<< "Min/max T:" << min(thermoScal.T()).value() << ' '
                         << max(thermoScal.T()).value() << endl;

Info << "El campo hscal calculado - max: " << max(hscal.primitiveField()) << "- min: " << min(hscal.primitiveField()) << endl;

// ************************************************************************* //

 //- fvm::laplacian((thermoScal.kappa()/(thermoScal.Cp() * thermoScal.rho())), hscal, "laplacian(kapparhoCp,T)")

