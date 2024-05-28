/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "Unal.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "phaseSystem.H"





#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"

#include "saturationModel.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureDiameterModels
{
    defineTypeNameAndDebug(Unal, 0);
    addToRunTimeSelectionTable
    (
        departureDiameterModel,
        Unal,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::
Unal::Unal
(
    const dictionary& dict
)
:
    departureDiameterModel(),
    phi_(readScalar(dict.lookup("phi")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureDiameterModels::
Unal::~Unal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureDiameterModels::
Unal::dDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{

    scalarField TsubUnal (Tsatw - Tl);
    scalarField TsupUnal (Tw - Tsatw);

    scalarField A((TsupUnal/(2.0*rhoLiquidw*L)) * pow(kw*rhoLiquidw*Cpw/pi,0.5));

    scalarField Mg (TsubUnal/(2.0*(1.0-rhoVaporw/rhoLiquidw)));

    scalarField bUnal =  Foam::neg(TsubUnal-3.0)*(Mg *exp(TsubUnal/3.0-1))
                       + Foam::pos(TsubUnal-3.0)*Mg; //Foam::neg0()


    //dimensionedScalar UbUnal("UbUnal", dimensionSet(0,1,-1,0,0,0,0), 0.61);

    //dDep_ = 2.42e-5 * pow(pw,0.709)*A;

    scalarField UwUnal (mag(Uw.patchInternalField()));


    scalarField phiUnal (max(pow(UwUnal/0.61,0.47)),1.0);

    scalarField DUnal (2.42e-5 * pow(pw,0.709)*A/(bUnal,pow(phiUnal,0.5)));



    return
         2.42e-5 * pow(p,0.709)* theta1 / sqrt(theta2*theta3);
        //0.0012*pow(rhoM, 0.9)*0.0208*phi_
       //*sqrt(sigmaw/(mag(g.value())*(rhoLiquid - rhoVapor)));
}


void Foam::wallBoilingModels::departureDiameterModels::
Unal::write(Ostream& os) const
{
    departureDiameterModel::write(os);
    writeEntry(os, "phi", phi_);
}


// ************************************************************************* //
