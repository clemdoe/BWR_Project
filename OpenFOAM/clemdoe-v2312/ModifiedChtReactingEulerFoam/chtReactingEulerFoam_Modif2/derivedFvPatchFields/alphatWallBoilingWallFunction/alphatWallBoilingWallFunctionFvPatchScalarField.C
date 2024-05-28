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

#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"
#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "saturationModel.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarField::phaseType,
    2
>::names[] =
{
    "vapor",
    "liquid"
};

const Foam::NamedEnum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarField::phaseType,
    2
>
Foam::compressible::
alphatWallBoilingWallFunctionFvPatchScalarField::phaseTypeNames_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    otherPhaseName_("vapor"),
    phaseType_(liquidPhase),
    relax_(0.5),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr)
{
    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    otherPhaseName_(dict.lookup("otherPhase")),
    phaseType_(phaseTypeNames_.read(dict.lookup("phaseType"))),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.5)),
    AbyV_(p.size(), 0),
    alphatConv_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    qq_(p.size(), 0),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr)
{

    // Check that otherPhaseName != this phase
    if (internalField().group() == otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapor phase that "
            << "corresponds to the liquid base of vice versa" << nl
            << "This phase: " << internalField().group() << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }

    switch (phaseType_)
    {
        case vaporPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            dmdt_ = 0;

            break;
        }
        case liquidPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            nucleationSiteModel_ =
                wallBoilingModels::nucleationSiteModel::New
                (
                    dict.subDict("nucleationSiteModel")
                );

            departureDiamModel_ =
                wallBoilingModels::departureDiameterModel::New
                (
                    dict.subDict("departureDiamModel")
                );

            departureFreqModel_ =
                wallBoilingModels::departureFrequencyModel::New
                (
                    dict.subDict("departureFreqModel")
                );

            if (dict.found("dDep"))
            {
                dDep_ = scalarField("dDep", dict, p.size());
            }

            if (dict.found("qQuenching"))
            {
                qq_ = scalarField("qQuenching", dict, p.size());
            }

            break;
        }
    }

    if (dict.found("alphatConv"))
    {
        alphatConv_ = scalarField("alphatConv", dict, p.size());
    }

    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(mapper(psf.alphatConv_)),
    dDep_(mapper(psf.dDep_)),
    qq_(mapper(psf.qq_)),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


alphatWallBoilingWallFunctionFvPatchScalarField::
alphatWallBoilingWallFunctionFvPatchScalarField
(
    const alphatWallBoilingWallFunctionFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_),
    AbyV_(psf.AbyV_),
    alphatConv_(psf.alphatConv_),
    dDep_(psf.dDep_),
    qq_(psf.qq_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    departureFreqModel_(psf.departureFreqModel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatWallBoilingWallFunctionFvPatchScalarField::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(otherPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
dmdt(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdt_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdt requested for invalid phasePair!"
            << abort(FatalError);

        return dmdt_;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarField::
mDotL(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return mDotL_;
    }
    else
    {
        FatalErrorInFunction
            << " mDotL requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}

void alphatWallBoilingWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Check that partitioningModel has been constructed
    if (!partitioningModel_.valid())
    {
        FatalErrorInFunction
            << "partitioningModel has not been constructed!"
            << abort(FatalError);
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const saturationModel& satModel =
        db().lookupObject<saturationModel>("saturationModel");

    const label patchi = patch().index();

//    Foam::Info << "Variable FLUID" << fluid << endl;
//    Foam::Info << "Variable satModel" << satModel.name() << endl;
//    Foam::Info << "Variable patchi" << patch().name() << endl;

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor
            (
                fluid.phases()[internalField().group()]
            );

            // Vapor Liquid phase fraction at the wall
            const scalarField vaporw(vapor.boundaryField()[patchi]);

            // NOTE! Assumes 1-thisPhase for liquid fraction in
            // multiphase simulations
            const scalarField fLiquid
            (
                partitioningModel_->fLiquid(1-vaporw)
            );

            operator==
            (
                calcAlphat(*this)*(1 - fLiquid)/max(vaporw, scalar(1e-8))
            );
            break;
        }
        case liquidPhase:
        {
            // Check that nucleationSiteModel has been constructed
            if (!nucleationSiteModel_.valid())
            {
                FatalErrorInFunction
                    << "nucleationSiteModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that departureDiameterModel has been constructed
            if (!departureDiamModel_.valid())
            {
                FatalErrorInFunction
                    << "departureDiameterModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that nucleationSiteModel has been constructed
            if (!departureFreqModel_.valid())
            {
                FatalErrorInFunction
                    << "departureFrequencyModel has not been constructed!"
                    << abort(FatalError);
            }

            const phaseModel& liquid
            (
                fluid.phases()[internalField().group()]
            );

            const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

            // Retrieve turbulence properties from models
            const phaseCompressibleTurbulenceModel& turbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                );
            const phaseCompressibleTurbulenceModel& vaporTurbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        vapor.name()
                    )
                );

            const nutWallFunctionFvPatchScalarField& nutw =
                nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

            const scalar Cmu25(pow025(nutw.Cmu()));

            const scalarField& y = turbModel.y()[patchi];

            const tmp<scalarField> tmuw = turbModel.mu(patchi);
            const scalarField& muw = tmuw();

            const tmp<scalarField> talphaw = liquid.thermo().alpha(patchi);
            const scalarField& alphaw = talphaw();

            const tmp<volScalarField> tk = turbModel.k();
            const volScalarField& k = tk();

            //Info << "Unidades k: " << k.dimensions() << endl;


            const fvPatchScalarField& kw = k.boundaryField()[patchi];

            const fvPatchVectorField& Uw =
                turbModel.U().boundaryField()[patchi];
            const scalarField magUp(mag(Uw.patchInternalField() - Uw));
            const scalarField magGradUw(mag(Uw.snGrad()));

            const fvPatchScalarField& rhoLiquidw =
                turbModel.rho().boundaryField()[patchi];

            const fvPatchScalarField& rhoVaporw =
                vaporTurbModel.rho().boundaryField()[patchi];

            const fvPatchScalarField& pw =
                liquid.thermo().p().boundaryField()[patchi];

            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];


            const scalarField Tc(Tw.patchInternalField());

           // Foam::Info << "El valor de Tc es: " << Tc << endl;

            const scalarField uTau(Cmu25*sqrt(kw));

            const scalarField yPlus(uTau*y/(muw/rhoLiquidw));


            //------------------------------------------//
            //Ojoo!!!! Aca mande mano al yPlus


//            scalarField kwmodif (kw*1e-8+average(kw));
//            const scalarField uTauModif(Cmu25*sqrt(kwmodif));
//            scalarField yPlusPesado (uTauModif*y/(muw/rhoLiquidw));
//            //------------------------------------------//


            // Aca solo voy a reportar el valor del yPlus para anotarlo
            //Foam::Info << "yPlus->" << "max: " << max(yPlus) << " min: "  << min(yPlus) << " Promedio: " << average(yPlus) << endl;


            const scalarField Pr(muw/alphaw);

            // Molecular-to-turbulent Prandtl number ratio
            const scalarField Prat(Pr/Prt_);

            // Thermal sublayer thicknessUnable to find turbulence model
            const scalarField P(this->Psmooth(Prat));

            const scalarField yPlusTherm(this->yPlusTherm(nutw, P, Prat));

            tmp<volScalarField> tCp = liquid.thermo().Cp();
            const volScalarField& Cp = tCp();
            const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

            // Saturation temperature
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());

            const volScalarField& Tsat = tTsat();
            const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
            const scalarField Tsatc(Tsatw.patchInternalField());

            const fvPatchScalarField& hew
                = liquid.thermo().he().boundaryField()[patchi];

            const scalarField hw
            (
                liquid.thermo().he().member() == "e"
              ? hew.patchInternalField() + pw/rhoLiquidw.patchInternalField()
              : hew.patchInternalField()
            );

            const scalarField L
            (
                vapor.thermo().he().member() == "e"
              ? vapor.thermo().he(pw, Tsatc, patchi) + pw/rhoVaporw - hw
              : vapor.thermo().he(pw, Tsatc, patchi) - hw
            );

            // Liquid phase fraction at the wall
            const scalarField liquidw(liquid.boundaryField()[patchi]);

            const scalarField fLiquid(partitioningModel_->fLiquid(liquidw));

            // Convective thermal diffusivity
            alphatConv_ = calcAlphat(alphatConv_);

            for (label i=0; i<10; i++)
            {
                // Liquid temperature at y+=250 is estimated from logarithmic
                // thermal wall function (Koncar, Krepper & Egorov, 2005)
                const scalarField Tplus_y250
                (
                    Prt_*(log(nutw.E()*250)/nutw.kappa() + P)
                );

                const scalarField Tplus
                (
                    Prt_*(log(nutw.E()*yPlus)/nutw.kappa() + P)
                );

               /* const scalarField TplusModif
                (
                    Prt_*(log(nutw.E()*150)/nutw.kappa() + P)
                );*/


                const scalarField Tl
                (
                   // Version original
                    max
                    (
                        Tc - 40,
                        Tw - (Tplus_y250/Tplus)*(Tw - Tc)
                    )

                   //Version modificada
                         /*   max
                            (
                                Tc - 40,
                                Tw - (max(Tplus_y250)/TplusModif)*(Tw - Tc)
                            )*/
                );



              //  Info << "Tl  - max: " << max(Tl) << " min: " << min(Tl) << endl;

               /* const scalarField Tlmodif
                       (
                            max
                            (
                                Tc - 40,
                                Tw - (max(Tplus_y250)/TplusModif)*(Tw - Tc)
                            )

                            );*/

//                Info << "Diferentes parametros"<<endl;
//                Info << "z[m]" <<"\t" << "Tl" << "\t" << "Tlmodif" << endl;
//                forAll(Tl,patchi)
//                {
//                 Info << patch().Cf()[patchi][2] << "\t"
//                      << yPlus[patchi] << "\t"
//                      << Tl[patchi] << "\t"
//                      << Tlmodif[patchi] << "\t"
//                      << Tplus[patchi] << "\t"
//                      << Tplus_y250[patchi] << endl;
//                }

                // Nucleation site density:
                const scalarField N
                (
                    nucleationSiteModel_->N
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl,
                        Tsatw,
                        L
                    )
                );

                // Bubble departure diameter:
                dDep_ = departureDiamModel_->dDeparture
                (
                    liquid,
                    vapor,
                    patchi,
                    Tl,
                    Tsatw,
                    L
                );

//                // -------------------------------------------------------------------- //
//                //Modificacion - Incorporacion Modelo de Unal

//                 const scalarField TsubUnal (Tsatw - Tl);
//                 const scalarField TsupUnal (Tw - Tsatw);


//////                //Caso A
//////                const scalarField TsubUnal (1.6+ Tsatw - Tl);
//////                const scalarField TsupUnal (4.7+ Tw - Tsatw);

////////                const scalarField TsubUnal (1.0+ Tsatw - Tl);
////////                const scalarField TsupUnal (2.5+ Tw - Tsatw);


//////                //Caso C
////////                const scalarField TsubUnal (2.2+ Tsatw - Tl);
////////                const scalarField TsupUnal (10.0+ Tw - Tsatw);



////                   // La conductividad del solido cobre es: 400
////                   // El rho del cobre es: 8960
////                   // El Cp cobre: 385

//                   //const scalarField A((TsupUnal/(2.0*rhoVaporw*L)) * pow(kw*rhoLiquidw*Cpw/pi,0.5));
//                   //const scalarField A((TsupUnal/(2.0*13.01*84730.0)) * pow(400.0*8960.0*385.0/pi,0.5));
//                   const scalarField A((TsupUnal/(2.0*13.01*L)) * pow(400.0*8960.0*385.0/pi,0.5));

//                   //const scalarField Mg (TsubUnal/(2.0*(1.0-rhoVaporw/rhoLiquidw)));
//                   const scalarField Mg (TsubUnal/(2.0*(1.0-13.01/rhoLiquidw)));

//                   const scalarField bUnal =  Foam::neg(TsubUnal-3.0)*(Mg *exp(TsubUnal/3.0-1.0))
//                                            + Foam::pos(TsubUnal-3.0)*Mg;

//                   const scalarField UwUnal (mag(Uw.patchInternalField()));

//                   const scalarField phiUnal (max(pow(UwUnal/0.61,0.47),1.0));

//                   const scalarField pwUnal((pw.patchInternalField())*0.0+101325);
//                   //const scalarField pwUnal (min(max((pw.patchInternalField()),pw.patchInternalField()*0.0+2.5e4),pw.patchInternalField()*0.0+1.2e5));

//                   //scalarField DUnal ((max(1.0e-6,(0.0000242 * pow(pwUnal,0.709)*A/(bUnal*pow(phiUnal,0.5))))));
//                   //const scalarField DUnal (min((max(1.0e-6,(0.0005+0.0000242 * pow(pwUnal,0.709)*A/(bUnal*pow(phiUnal,0.5))))),0.0012));
//                   const scalarField DUnal (min((max(1.0e-6,(0.0000242 * pow(pwUnal,0.709)*A/(bUnal*pow(phiUnal,0.5))))),0.003));

//                   dDep_ = 1.0*DUnal;

////                   Info << "Reporto los flujos de calor"<<endl;
////                                  Info << "z[m]" <<"\t" << "TsubUnal" << " \t \t" << "TsupUnal" << "\t" << "A" << "\t" << "bUnal" << endl;
////                                  forAll(TsubUnal,patchi)
////                                  {
////                                   Info << patch().Cf()[patchi][2] << "\t \t"
////                                        << TsubUnal[patchi] << "\t"
////                                        << TsupUnal[patchi]  << "\t"
////                                        << A[patchi]          << "\t"
////                                        << UwUnal[patchi]          <<  "\t"
////                                        << bUnal[patchi]          <<  "\t"
////                                        << phiUnal[patchi]         <<  "\t"
////                                        << dDep_[patchi]          << "\t"
////                                        << pwUnal[patchi] <<  endl;

////                                  }

               //     Info << "DUnal - Max: " << max(dDep_) << " - Min: " <<  min(dDep_) << " - Average: " << average(dDep_) << endl;
//                    Info << "pUnal - Max: " << max(pwUnal) << " - Min: " <<  min(pwUnal) << " - Average: " << average(pwUnal) << endl;


                // -------------------------------------------------------------------- //


                // Bubble departure frequency:
                const scalarField fDep
                (
                    departureFreqModel_->fDeparture
                    (
                        liquid,
                        vapor,
                        patchi,
                        dDep_
                    )
                );

                // Area fractions:

                // Del Valle & Kenning (1985)
                const scalarField Ja
                (
                    rhoLiquidw*Cpw*(Tsatw - Tl)/(rhoVaporw*L)
                );

                const scalarField Al
                (
                    fLiquid*4.8*exp( min(-Ja/80, log(vGreat)))
                );

                const scalarField A2(min(pi*sqr(dDep_)*N*Al/4, scalar(1)));
                const scalarField A1(max(1 - A2, scalar(1e-4)));
                const scalarField A2E(min(pi*sqr(dDep_)*N*Al/4, scalar(5)));
                //El A2E que vale es el de arriba, esto es una prueba
                //const scalarField A2E(pi*sqr(dDep_)*N*Al/4);


                //Peso la areas por un factor scal
                //scalar scal (0.2);
//                const scalarField A2(min(pi*sqr(dDep_)*N*Al/4, scalar(1))*(1+scal+0.1));
//                const scalarField A1(max(1 - A2, scalar(1e-4))*(1-scal));
//                const scalarField A2E(min(pi*sqr(dDep_)*N*Al/4, scalar(5))*(1+0.2));


                // Volumetric mass source in the near wall cell due to the
                // wall boiling

                dmdt_ =
                    (1 - relax_)*dmdt_
                  + relax_*(1.0/6.0)*A2E*dDep_*rhoVaporw*fDep*AbyV_;

                // Volumetric source in the near wall cell due to the wall
                // boiling
                mDotL_ = dmdt_*L;

                // Quenching heat transfer coefficient
                const scalarField hQ
                (
                    2*(alphaw*Cpw)*fDep*sqrt((0.8/fDep)/(pi*alphaw/rhoLiquidw))

                );

//                const scalarField hQmodif
//                (
//                    2*alphaw / sqrt((pi/fDep)*(alphaw/(rhoLiquidw*Cpw)))

//                );

//                forAll(hQ,patchi)
//                {
//                Info << patch().Cf()[patchi][2] << "\t"
//                     << hQ[patchi] << "\t"
//                     << hQmodif[patchi] << endl;
//                }


                // Quenching heat flux
                qq_ = (A2*hQ*max(Tw - Tl, scalar(0)));

                // Effective thermal diffusivity that corresponds to the
                // calculated convective, quenching and evaporative heat fluxes

                operator==
                (
                    (
                        A1*alphatConv_
                      + (qq_ + qe())/max(hew.snGrad(), scalar(1e-16))
                    )
                   /max(liquidw, scalar(1e-8))
                );

                scalarField TsupPrev(max((Tw - Tsatw), scalar(0)));
                const_cast<fvPatchScalarField&>(Tw).evaluate();
                scalarField TsupNew(max((Tw - Tsatw), scalar(0)));

                scalar maxErr(max(mag(TsupPrev - TsupNew)));

//                Info << "Tw " << Tw << endl;
//                Info << "Tsatw " << Tsatw << endl;

                // ------------------------------------------------------------- //
                    // Ojo porque aca paso todo a W/cm2 que son como los reporta el paper
                const scalarField qee
                        (
                           (fLiquid*qe())/10000
                            //(pi/6) * (dDep_)*(dDep_)*(dDep_) * fDep * N * rhoVaporw *84730/10000
                            //(pi/6) * (dDep_)*(dDep_)*(dDep_) * fDep * N * rhoVaporw * L/10000
                         );

                const scalarField qc
                (
                    fLiquid*A1*(alphatConv_ + alphaw)*hew.snGrad()/10000
                );

                const scalarField qEff
                (
                    liquidw*(*this + alphaw)*hew.snGrad()/10000
                );

                const scalarField qquen
                (
                  (fLiquid*qq())/10000
                );



                Info << "Reporto los flujos de calor"<<endl;
                Info << "z[m]" <<"\t" << "Qtotal" << " \t \t" << "Qquenching" << "\t" << "Qconvection" << "\t" << "Qevapor" << endl;
                forAll(qEff,patchi)
                {
                 Info << patch().Cf()[patchi][2] << "\t"
                      << qEff[patchi] << "\t"
                      << qquen[patchi]<< "\t"
                      << qc[patchi]   << "\t"
                      << qee[patchi]   << "\t"
                      << dDep_[patchi]   << "\t"
                      << Tw[patchi]   <<  "\t"
                      << N[patchi]   <<  "\t"
                      << fDep[patchi]   <<  "\t"
                      << Tl[patchi] << "\t"
                      << rhoVaporw[patchi] << "\t"
                         << fLiquid[patchi] <<   endl;
                      //<< DUnal[patchi]<<   endl;
                }

                // ------------------------------------------------------------- //

                Info << "Num. Iter bucle Newton: " << i << " - El error: " << maxErr << endl;


                if (debug)
                {
                    const scalarField qc
                    (
                        fLiquid*A1*(alphatConv_ + alphaw)*hew.snGrad()
                    );

                    const scalarField qEff
                    (
                        liquidw*(*this + alphaw)*hew.snGrad()
                    );

                    Info<< "  L: " << gMin(L) << " - " << gMax(L) << endl;
                    Info<< "  Tl: " << gMin(Tl) << " - " << gMax(Tl) << endl;
                    Info<< "  N: " << gMin(N) << " - " << gMax(N) << endl;
                    Info<< "  dDep_: " << gMin(dDep_) << " - "
                        << gMax(dDep_) << endl;
                    Info<< "  fDep: " << gMin(fDep) << " - "
                        << gMax(fDep) << endl;
                    Info<< "  Al: " << gMin(Al) << " - " << gMax(Al) << endl;
                    Info<< "  A1: " << gMin(A1) << " - " << gMax(A1) << endl;
                    Info<< "  A2: " << gMin(A2) << " - " << gMax(A2) << endl;
                    Info<< "  A2E: " << gMin(A2E) << " - "
                        << gMax(A2E) << endl;
                    Info<< "  dmdtW: " << gMin(dmdt_) << " - "
                        << gMax(dmdt_) << endl;
                    Info<< "  qc: " << gMin(qc) << " - " << gMax(qc) << endl;
                    Info<< "  qq: " << gMin(fLiquid*qq()) << " - "
                        << gMax(fLiquid*qq()) << endl;
                    Info<< "  qe: " << gMin(fLiquid*qe()) << " - "
                        << gMax(fLiquid*qe()) << endl;
                    Info<< "  qEff: " << gMin(qEff) << " - "
                        << gMax(qEff) << endl;
                    Info<< "  alphat: " << gMin(*this) << " - "
                        << gMax(*this) << endl;
                    Info<< "  alphatConv: " << gMin(alphatConv_)
                        << " - " << gMax(alphatConv_) << endl;
                }

                if (maxErr < 1e-1)
                {
                    if (i > 0)
                    {
                        Info<< "Wall boiling wall function iterations: "
                            << i + 1 << endl;
                    }
                    break;
                }

            }
            break;



        }
        default:
        {
            FatalErrorInFunction
                << "Unknown phase type. Valid types are: "
                << phaseTypeNames_ << nl << exit(FatalError);
        }


    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    writeEntry(os, "phaseType", phaseTypeNames_[phaseType_]);

    writeEntry(os, "relax", relax_);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;
            break;
        }
        case liquidPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("nucleationSiteModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            nucleationSiteModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureDiamModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureDiamModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureFreqModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureFreqModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            break;
        }
    }

    writeEntry(os, "otherPhase", otherPhaseName_);
    writeEntry(os, "dmdt", dmdt_);
    writeEntry(os, "dDep", dDep_);
    writeEntry(os, "qQuenching", qq_);
    writeEntry(os, "alphatConv", alphatConv_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatWallBoilingWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
