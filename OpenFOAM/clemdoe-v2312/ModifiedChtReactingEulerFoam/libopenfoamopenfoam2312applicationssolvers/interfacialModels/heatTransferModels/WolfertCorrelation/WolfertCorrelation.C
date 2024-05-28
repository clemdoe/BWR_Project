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

\*---------------------------------------------------------------------------*/

#include "WolfertCorrelation.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(WolfertCorrelation, 0);
    addToRunTimeSelectionTable(heatTransferModel, WolfertCorrelation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::WolfertCorrelation::WolfertCorrelation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair),
    Cpfase_("Cpfase", dimensionSet(0, 2, -2, -1, 0,0,0), dict),
    Cc_(dict.lookupOrDefault<scalar>("Cc", 2.0))

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::WolfertCorrelation::~WolfertCorrelation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::WolfertCorrelation::K(const scalar residualAlpha) const
{

    // Cpfase es el Cp de la fase contina

    volScalarField magUr  (mag(pair_.Ur()));

    volScalarField kTurbcont (pair_.continuous().kappaEff() - pair_.continuous().kappa());

    volScalarField kcont (pair_.continuous().kappa());

    volScalarField rhocon (pair_.continuous().rho());

    //Le puse un factor de escala 2

    volScalarField hpel (Cc_*rhocon * Cpfase_ *  sqrt( (3.1415926/4.0)*
                                                   (magUr/pair_.dispersed().d())
                                                   * (kcont/(rhocon * Cpfase_))
                                                   * (1.0/(kTurbcont/kcont +1.0))));

    return
        (6.0 * hpel *max(pair_.dispersed(), residualAlpha))/ pair_.dispersed().d();

}


// ************************************************************************* //
