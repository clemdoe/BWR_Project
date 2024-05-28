/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"


// ------ Esto lo incluyo para leer las propiedades de los dos fluidos ---- //
//#include "fixedMultiPhaseHeatFluxFvPatchScalarField.H"
#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase_Modif2(patch(), "undefined", "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    TwoFluidAcoupled_(false),
    TypeRegion_("unidefined-TypeRegion")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase_Modif2(patch(), ptf),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
    contactRes_(ptf.contactRes_),
    TwoFluidAcoupled_(ptf.TwoFluidAcoupled_),
    TypeRegion_(ptf.TypeRegion_)
{}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase_Modif2(patch(), dict),
    TnbrName_(dict.lookup("Tnbr")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),
    TwoFluidAcoupled_(readBool(dict.lookup("TwoFluidAcoupled"))),
    TypeRegion_(dict.lookup("TypeRegion"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::
turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
(
    const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase_Modif2(patch(), wtcsf),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    kappaLayers_(wtcsf.kappaLayers_),
    contactRes_(wtcsf.contactRes_),
    TwoFluidAcoupled_(wtcsf.TwoFluidAcoupled_),
    TypeRegion_(wtcsf.TypeRegion_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    // Calculate the temperature by harmonic averaging
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    typedef turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2 thisType;

    //Si es cero resuelve una acoplamiento fluido/fluido o fluido/solido
    //Esto se define en la condicion de borde dentro del campo T
    if( TwoFluidAcoupled_ == 0)
    {
        Info << "Region: " << TypeRegion_ << endl;

    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }

    const thisType& nbrField = refCast<const thisType>(nbrTp);

    // Swap to obtain full local values of neighbour internal field
    tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
    tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

    if (contactRes_ == 0.0)
    {
        nbrIntFld.ref() = nbrField.patchInternalField();
        nbrKDelta.ref() = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        nbrIntFld.ref() = nbrField;
        nbrKDelta.ref() = contactRes_;
    }

    mpp.distribute(nbrIntFld.ref());
    mpp.distribute(nbrKDelta.ref());

    tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();

//    Info << "------- Calculando la Armonica de la Temperatura ---------- " << endl;
//    Info << "El valor de kappa en esta region: " << kappa(*this) << endl;
//    Info << "El valor de deltaCoeffs en esta region: " << patch().deltaCoeffs() << endl;
//    Info << "El valor de myKDelta: " << kappa(*this)*patch().deltaCoeffs() << endl;
//    Info << "El valor de nbrKDelta: " << nbrKDelta() << endl;
//    Info << "El valor de mixFraction: " << nbrKDelta()/(nbrKDelta() + myKDelta()) << endl;
//    Info << "El valor de this " << *this << endl;
//    Info << "---------------------------------------------------------" << endl;


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

//    Info << "Region Fluido/Solido o Solido/Fluido " << endl;
//    Info << "nbrKDelta(): "
//         << " - max: " <<   max(nbrKDelta())
//         << " - min: " <<   min(nbrKDelta()) << endl;
//    Info << "myKDelta(): "
//         << " - max:" << max(myKDelta())
//         << " - min:" << min(myKDelta())<< endl;
//    Info << "valueFraction(): "
//         << " - max: " << max(nbrKDelta()/(nbrKDelta() + myKDelta()))
//         << " - min: " << min(nbrKDelta()/(nbrKDelta() + myKDelta()))<< endl;


    this->refValue() = nbrIntFld();
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

    mixedFvPatchScalarField::updateCoeffs();

    //----------------------------------------------------//
    // Ojo porque esto lo saque de la parte que decia debug.
    // despues lo tengo q comentar
    scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

    Info << "Patch name: " << patch().name() << " -- Q [W]: " << Q << endl;



    //Esto lo use para reportar W/cm2 pero sobre todo el pacht, y no como una suma
    //scalarField Qregion = kappa(*this)*snGrad()/10000;
    //Info << "Patch name: " << patch().name() << " -- Q [W]: " << Qregion << endl;


    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    }

    // ---------------------------------------------------------------------------------------- //

    //Si es TwoFluidAcoupled es 1, entonces resuelve el acoplamietno entre una region monofasica con una de dos fluidos
    //Esto se define en la condicion de borde dentro del campo T
    else if(TwoFluidAcoupled_ == 1)
    {
        if(TypeRegion_ == "TwoFluid")
        {
         //   Info << "Region: " << TypeRegion_ << endl;

            const fvPatchScalarField& nbrTp =
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

            if (!isA<thisType>(nbrTp))
            {
                FatalErrorInFunction
                    << "Patch field for " << internalField().name() << " on "
                    << patch().name() << " is of type " << thisType::typeName
                    << endl << "The neighbouring patch field " << TnbrName_ << " on "
                    << nbrPatch.name() << " is required to be the same, but is "
                    << "currently of type " << nbrTp.type() << exit(FatalError);
            }

            const thisType& nbrField = refCast<const thisType>(nbrTp);

            // Swap to obtain full local values of neighbour internal field
            tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
            tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

            if (contactRes_ == 0.0)
            {
                nbrIntFld.ref() = nbrField.patchInternalField();
                nbrKDelta.ref() = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
            }
            else
            {
                nbrIntFld.ref() = nbrField;
                nbrKDelta.ref() = contactRes_;
            }

            mpp.distribute(nbrIntFld.ref());
            mpp.distribute(nbrKDelta.ref());

            //En esta parte calculo alpha_vapor * kappaEff_vapor * DeltaCoeff + alpha_liq * kappaEff_liq * DeltaCoeff

            // Lookup the fluid model
            const phaseSystem& fluid =
                (
                    db().lookupObject<phaseSystem>("phaseProperties")
                );

            const scalarField& Tp = *this;

            scalarField B(Tp.size(), scalar(0));

            forAll(fluid.phases(), phasei)
            {
                const phaseModel& phase = fluid.phases()[phasei];

                const fvPatchScalarField& alpha =
                    phase.boundaryField()[patch().index()];

                const scalarField kappaEff(phase.kappaEff(patch().index()));

                B += alpha*kappaEff*patch().deltaCoeffs();

            }

            tmp<scalarField> myKDelta =  B*1.0;


//            Info << "Region Dos-Fluidos/Solido " << endl;
//            Info << "nbrKDelta(): "
//                 << " - max: " <<   max(nbrKDelta())
//                 << " - min: " <<   min(nbrKDelta()) << endl;
//            Info << "myKDelta(): "
//                 << " - max:" << max(myKDelta())
//                 << " - min:" << min(myKDelta())<< endl;
//            Info << "valueFraction(): "
//                 << " - max: " << max(nbrKDelta()/(nbrKDelta() + myKDelta()))
//                 << " - min: " << min(nbrKDelta()/(nbrKDelta() + myKDelta()))<< endl;

            this->refValue() = nbrIntFld();
            this->refGrad() = 0.0;
            this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

            mixedFvPatchScalarField::updateCoeffs();

            //----------------------------------------------------//
            scalar Q = gSum((myKDelta/patch().deltaCoeffs())*patch().magSf()*snGrad());

        //    Info << "Patch name: " << patch().name() << " -- Q [W]: " << Q << endl;

            // Restore tag
            UPstream::msgType() = oldTag;
        }
        else if(TypeRegion_ == "solido")
        {

            const fvPatchScalarField& nbrTp =
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

           // Info << "Region: " << TypeRegion_ << endl;

            const thisType& nbrField = refCast<const thisType>(nbrTp);

//            // Swap to obtain full local values of neighbour internal field
            tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
            tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

            if (contactRes_ == 0.0)
            {
                nbrIntFld.ref() = nbrField.patchInternalField();
                nbrKDelta.ref() = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
            }
            else
            {
                nbrIntFld.ref() = nbrField;
                nbrKDelta.ref() = contactRes_;
            }

            mpp.distribute(nbrIntFld.ref());
            mpp.distribute(nbrKDelta.ref());

            tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();

//            Info << "Region Solido/Dos-Fluidos " << endl;
//            Info << "nbrKDelta(): "
//                 << " - max: " <<   max(nbrKDelta())
//                 << " - min: " <<   min(nbrKDelta()) << endl;
//            Info << "myKDelta(): "
//                 << " - max:" << max(myKDelta())
//                 << " - min:" << min(myKDelta())<< endl;
//            Info << "valueFraction(): "
//                 << " - max: " << max(nbrKDelta()/(nbrKDelta() + myKDelta()))
//                 << " - min: " << min(nbrKDelta()/(nbrKDelta() + myKDelta()))<< endl;


            this->refValue() = nbrIntFld();
            this->refGrad() = 0.0;
            this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

            mixedFvPatchScalarField::updateCoeffs();

            //----------------------------------------------------//
            scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

            //Info << "Patch name: " << patch().name() << " -- Q[W]: " << Q << endl;

            // Restore tag
            UPstream::msgType() = oldTag;
        }
        else
        {

        FatalErrorInFunction
            << "No esta implementado obtener la armonica de la temperatura para el modelo fluidos monofasico con dos-fluidos"
            << exit(FatalError);
        }

    }
    // ---------------------------------------------------------------------------------------- //



}


void turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);
    writeEntry(os, "TwoFluidAcoupled", TwoFluidAcoupled_);
    writeEntry(os, "TypeRegion", TypeRegion_);

    temperatureCoupledBase_Modif2::write(os);
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
