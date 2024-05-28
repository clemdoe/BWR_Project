/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "hyperbolicTsubDiameter.H"
#include "phaseSystem.H"
#include "saturationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(hyperbolicTsub, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        hyperbolicTsub,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::hyperbolicTsub::hyperbolicTsub
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    liquidPhaseName_(diameterProperties.lookup("liquidPhase")),
    d2_("d2", dimLength, diameterProperties.lookupOrDefault("d2", 0.00015)),
    Tsub2_
    (
        "Tsub2",
         dimTemperature,
         diameterProperties.lookupOrDefault("Tsub2", -5.0)
    ),
    d1_
    (
        "d1",
        dimLength,
        diameterProperties.lookupOrDefault("d1", 0.002)
    ),
    Tsub1_
    (
        "Tsub1",
        dimTemperature,
        diameterProperties.lookupOrDefault("Tsub1", 13.5)
    ),

    TsubScale_
    (
        "TsubScale",
        dimTemperature,
        diameterProperties.lookupOrDefault("TsubScale", 4.0)
    ),

    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh(),
        d1_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::hyperbolicTsub::~hyperbolicTsub()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::hyperbolicTsub::d() const
{
    return d_;
}


void Foam::diameterModels::hyperbolicTsub::correct()
{
    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            phase_.mesh().lookupObject<phaseSystem>("phaseProperties")
        );

    const phaseModel& liquid(fluid.phases()[liquidPhaseName_]);

    if (phase_.mesh().foundObject<saturationModel>("saturationModel"))
    {
        const saturationModel& satModel =
            phase_.mesh().lookupObject<saturationModel>("saturationModel");


        //Asi esta implementado en ANSYS

        const volScalarField TsubL
        (
            satModel.Tsat(liquid.thermo().p()) - liquid.thermo().T()
        );



        const dimensionedScalar TsubMean
        (
            0.5 * (Tsub1_ + Tsub2_)
        );

        const volScalarField TanhArgum
        (
            (TsubL - TsubMean)/TsubScale_
        );

        const volScalarField TanhBlend
        (
             1.0 / (1.0 + exp(TanhArgum))
        );

        //d_ = TanhBlend * (d1_ - d2_)+d2_;

        /*d_ = max
        (
            d1_,
            min
            (
                d2_,
                (d1_*(Tsub - Tsub2_) + d2_*(Tsub - Tsub1_))/(Tsub2_ - Tsub1_)
            )
        );*/

        d_ = max
        (
            d1_,
            min
            (
                d2_,
               TanhBlend * (d1_ - d2_)+d2_
            )
        );


    }
}


bool Foam::diameterModels::hyperbolicTsub::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);
    diameterProperties_.lookup("liquidPhase") >> liquidPhaseName_;
    diameterProperties_.lookup("d2") >> d2_;
    diameterProperties_.lookup("Tsub2") >> Tsub2_;
    diameterProperties_.lookup("d1") >> d1_;
    diameterProperties_.lookup("Tsub1") >> Tsub1_;

    return true;
}


// ************************************************************************* //
