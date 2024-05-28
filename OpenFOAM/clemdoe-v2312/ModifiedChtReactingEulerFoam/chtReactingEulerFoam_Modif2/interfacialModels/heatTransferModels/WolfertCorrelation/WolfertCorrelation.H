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

Class
    Foam::heatTransferModels::WolfertCorrelation

Description
    Ranz-Marshall correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.

    The heat transfer coefficient is calculated using Wolfert’s correlation.
    Is used in the NEPTUNE code instead of the Wolfert’s correlation,
    has also been implemented but has not been used here. Then the condensation rate
    field is created from the heat transfer coefficient, cutting off potential unphysical negative values.


SourceFiles
    WolfertCorrelation.C

\*---------------------------------------------------------------------------*/

#ifndef WolfertCorrelation_H
#define WolfertCorrelation_H

#include "heatTransferModel.H"
#include "phaseSystem.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class WolfertCorrelation Declaration
\*---------------------------------------------------------------------------*/

class WolfertCorrelation
:
    public heatTransferModel
{
    // Private data

        //- Calor especifico de la fase
        const dimensionedScalar Cpfase_;

        //- Coefficiente de condensacion
            scalar Cc_;

public:

    //- Runtime type information
    TypeName("WolfertCorrelation");


    // Constructors

        //- Construct from components
        WolfertCorrelation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~WolfertCorrelation();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
