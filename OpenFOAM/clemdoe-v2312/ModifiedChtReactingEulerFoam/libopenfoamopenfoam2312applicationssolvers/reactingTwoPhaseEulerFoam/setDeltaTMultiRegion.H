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

Global
    setDeltaTMultiregion

Description
	Es igual al setDeltT pero se lee no solo la region de dos fluidos sino
	tambien la region de fluido monofasico y solido
	
    Restablezca el paso de tiempo para mantener un número de Courant máximo constante. 
	La reducción del intervalo de tiempo es inmediata, pero el aumento 
	se amortigua para evitar oscilaciones inestables.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
	//Region solido y fluido monofasico
    if (CoNumFluid == -great)
    {
        CoNumFluid = small;
    }

    if (DiNum == -great)
    {
        DiNum = small;
    }

    scalar maxDeltaTFluid = maxCo/(CoNumFluid + small);
    scalar maxDeltaTSolid = maxDi/(DiNum + small);

    scalar deltaTFluid =
        min
        (
            min(maxDeltaTFluid, 1.0 + 0.1*maxDeltaTFluid),
            1.2
        );

    //Esto es de la region dos-fluidos
 	scalar maxDeltaTFact = maxCo/(CoNum + small);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

	//Seteo el paso de tiempo tomando el minimo.
	runTime.setDeltaT
    (
        min
        (
            min(	
				min(deltaTFluid, maxDeltaTSolid),
					deltaTFact)*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //

