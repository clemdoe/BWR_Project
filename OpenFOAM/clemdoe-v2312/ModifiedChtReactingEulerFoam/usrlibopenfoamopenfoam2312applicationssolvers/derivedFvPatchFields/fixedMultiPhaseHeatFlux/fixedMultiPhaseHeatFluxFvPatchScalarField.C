/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "fixedMultiPhaseHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    relax_(1.0),
    Tmin_(0.0)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    Tmin_(dict.lookupOrDefault<scalar>("Tmin", 273))
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    q_(mapper(psf.q_)),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf
)
:
    fixedValueFvPatchScalarField(psf),
    q_(psf.q_),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}


Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::
fixedMultiPhaseHeatFluxFvPatchScalarField
(
    const fixedMultiPhaseHeatFluxFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    q_(psf.q_),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const scalarField& Tp = *this;


    scalarField A(Tp.size(), scalar(0));
    scalarField B(Tp.size(), scalar(0));
    scalarField Q(Tp.size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const fluidThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        const fvPatchScalarField& T =
            thermo.T().boundaryField()[patch().index()];

        const scalarField kappaEff(phase.kappaEff(patch().index()));

        if (debug)
        {
            scalarField q0(T.snGrad()*alpha*kappaEff);
            Q += q0;

            Info<< patch().name() << " " << phase.name()
                << ": Heat flux " << gMin(q0) << " - " << gMax(q0) << endl;
        }

        A += T.patchInternalField()*alpha*kappaEff*patch().deltaCoeffs();
        B += alpha*kappaEff*patch().deltaCoeffs();


        // ------------------------------------------------------ //

    //Busco las coordenadas de las caras de bordes

//        const vectorField Cb = patch().Cf();

//        const scalarField zCb = Cb.component(2);

//        Info<< "Field zCb:" << zCb << endl;

        //Info<< "Field cf: " << cf << endl;
        /*label celli(0);
        Info<< "All cells Cf():" << endl;
        forAll(Cb, facei) {
            Info<< cf[facei] << endl;
        }*/

        /*const vectorField& cf = patch().Cf();



        const scalarField& zcf = cf.component(2);
         Info<< "Field cf: " << cf << endl;
         Info<< "Field zcf:" << zcf << endl;

         label celli(0);
         const scalar z0 = cf[celli].z();
         Info<< "Just the first cell: " << z0 << endl;
         Info<< "Field z0: " << z0 << endl;

         Info<< "All cells Cf():" << endl;
         forAll(cf, facei) {
             Info<< cf[facei] << endl;
         }

         Info<< "Field zcf:" << zcf << endl;
*/
        // ------------------------------------------------------ //

    }

    if (debug)
    {
        Info<< patch().name() << " " << ": overall heat flux "
            << gMin(Q) << " - " << gMax(Q) << " W/m2, power: "
            << gSum(patch().magSf()*Q) << " W" << endl;
    }


    const vectorField Cb = patch().Cf();

    const scalarField zCb = Cb.component(2);

    //Info<< "Field zCb:" << zCb << endl;



    //Foam::Info << "El valor de Tp es: " << Tp  << endl;
   // Foam::Info << "El valor de q es: " << q_  << endl;
    //Foam::Info << "El valor de Tp es: " << Tp  << endl;
//    Foam::Info << "El valor de Tp es: " << Tp << T << endl;

    //Foam::Info << "El valor de q es: " << q_  << endl;

    //q_ = 1e-9;
    //q_ = 570000;

    //Potencia lineal maximo valor en z = 2 y valor cero en z = 0
    //q_ = q_ * (0.5*zCb);

    //Potencia lineal maximo valor en z = 0 y valor cero en z = 2
    //q_ = 570000;
    //q_ = q_ * (-zCb+2);
    //Foam::Info << "El valor de q modificado por z es: " << q_  << endl;

    //Porencia parabolica con maximo valor en z = 1 y valores ceros en z = 0 y z = 2
    //Potencia total equivalente al caso constante
    //q_ = 855000;
    //Info << "Aplico un potencia parabolica q_ * (-1.0*sqr(zCb) + 2.0 * zCb), con q_=855000 " << endl;
    //q_ = q_ * (-1.0* sqr(zCb) + 2.0 * zCb);
    //Foam::Info << "El valor de q modificado por z es: " << q_  << endl;

    //Foam::Info << "La posicion de los centros de cara: " << zCb << endl;




    //Foam::Info << "El valor de q_modif es: " << q_  << endl;

    //Esta es la version modificada sin relajacion
    //operator== (max(Tmin_,(q_ + A)/(B)));
    //Esta es la version original con relajacion
    operator==((1 - relax_)*Tp + relax_*max(Tmin_,(q_ + A)/(B)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fixedMultiPhaseHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "relax", relax_);
    writeEntry(os, "q", q_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedMultiPhaseHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
