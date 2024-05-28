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

#include "temperatureCoupledBase_Modif2.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "turbulentFluidThermoModel.H"


// ------ Esto lo incluyo para leer las propiedades de los dos fluidos ---- //
//#include "fixedMultiPhaseHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"



// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::temperatureCoupledBase_Modif2::KMethodType,
        5
    >::names[] =
    {
        "fluidThermo",
        "solidThermo",
        "directionalSolidThermo",
        "lookup",
        "TwoFluid"
    };
}


const Foam::NamedEnum<Foam::temperatureCoupledBase_Modif2::KMethodType, 5>
    Foam::temperatureCoupledBase_Modif2::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase_Modif2::temperatureCoupledBase_Modif2
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName)
{}


Foam::temperatureCoupledBase_Modif2::temperatureCoupledBase_Modif2
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("kappaMethod"))),
    kappaName_(dict.lookupOrDefault<word>("kappa", "none")),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAni","Anialpha"))
{
    switch (method_)
    {
        case mtDirectionalSolidThermo:
        {
            if (!dict.found("alphaAni"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'alphaAni'"
                       " required for 'kappaMethod' "
                    << KMethodTypeNames_[method_]
                    << exit(FatalIOError);
            }

            break;
        }

        case mtLookup:
        {
            if (!dict.found("kappa"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'kappa'"
                       " required for 'kappaMethod' "
                    <<  KMethodTypeNames_[method_] << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                       " or volSymmTensorField"
                    << exit(FatalIOError);
            }

            break;
        }

        default:
            break;
    }
}


Foam::temperatureCoupledBase_Modif2::temperatureCoupledBase_Modif2
(
    const fvPatch& patch,
    const temperatureCoupledBase_Modif2& base
)
:
    patch_(patch),
    method_(base.method_),
    kappaName_(base.kappaName_),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase_Modif2::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();


    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            word turbName(turbulenceModel::propertiesName);

//          Info << "turbName: " << turbName << endl;
//          Info << "basicThermo::dictName: " << basicThermo::dictName << endl;
//          Info << "mesh.name(): " << mesh.name() << endl;
//          Info << "mesh.nCells(): " << mesh.nCells() << endl;
//          Info << "Si existe objeto TurbulenceModel 1 - turbkappaEff: " << mesh.foundObject<turbulenceModel>(turbName) << endl;
//          Info << "Si existe objeto FluidoThermo 2 - thermo.kappa: " << mesh.foundObject<fluidThermo>(basicThermo::dictName) << endl;

            if
            (
                mesh.foundObject<turbulenceModel>(turbName)
            )
            {
                const turbulenceModel& turbModel =
                    mesh.lookupObject<turbulenceModel>(turbName);

//                Info << "kappa primario - max: " << max(turbModel.kappaEff(patchi))
//                     << " - min: " << min(turbModel.kappaEff(patchi)) << endl;

                return turbModel.kappaEff(patchi);
            }

            else if (mesh.foundObject<fluidThermo>(basicThermo::dictName))
            {
                const fluidThermo& thermo =
                    mesh.lookupObject<fluidThermo>(basicThermo::dictName);

                return thermo.kappa(patchi);
            }

            else
            {
                FatalErrorInFunction
                    << "kappaMethod defined to employ "
                    << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {

            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

//            Info << "kappa solid - max: " << max(thermo.kappa(patchi))
//                 << " - min: " << min(thermo.kappa(patchi)) << endl;

            return thermo.kappa(patchi);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const scalarField& pp = thermo.p().boundaryField()[patchi];

            const symmTensorField kappa(alphaAni*thermo.Cp(pp, Tp, patchi));

            const vectorField n(patch_.nf());

            return n & kappa & n;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(kappaName_))
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        kappaName_
                    );

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                       " or volSymmTensorField."
                    << exit(FatalError);
            }

            break;
        }

        //Caso seteado para el dos fluidos
        case TwoFluid:
        {
         // Lookup the fluid model
         const phaseSystem& fluid =
             (
                 mesh.lookupObject<phaseSystem>("phaseProperties")
             );

         scalarField KappaAlpha(patch_.size(),scalar(0));

         forAll(fluid.phases(), phasei)
         {
             const phaseModel& phase = fluid.phases()[phasei];

             const fvPatchScalarField& alphaRef =
                phase.boundaryField()[patch_.index()];

             const scalarField kappaEffRef(phase.kappaEff(patch_.index()));
             //Es lo mismo patchi = patch_.index()
             //const scalarField kappaEffRef(phase.kappaEff(patch_.index()));

//             Info << "Water and vapor " << endl;
//             Info << "alpha." << phase.name()
//                  << " - max: " << max(alphaRef)
//                  << " - min: " << min(alphaRef) << endl;
//             Info << "kappaEff." << phase.name()
//                  << " - max: " << max(kappaEffRef)
//                  << " - min: " << min(kappaEffRef) << endl;

             KappaAlpha += alphaRef*kappaEffRef*1.0;
             //KappaAlpha += kappaEffRef*1.0;
         }

//         Info << "kappa_v*alpha_v + kappa_l*alpha_l :"
//              << " max: " << max(KappaAlpha)
//              << " min: " << min(KappaAlpha) << endl;

         return KappaAlpha*1.0;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "    Please set 'kappaMethod' to one of "
                << KMethodTypeNames_.toc()
                << " and 'kappa' to the name of the volScalar"
                << " or volSymmTensor field (if kappa=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::temperatureCoupledBase_Modif2::write(Ostream& os) const
{
    writeEntry(os, "kappaMethod", KMethodTypeNames_[method_]);
    writeEntry(os, "kappa", kappaName_);
    writeEntry(os, "alphaAni", alphaAniName_);
}


// ************************************************************************* //
