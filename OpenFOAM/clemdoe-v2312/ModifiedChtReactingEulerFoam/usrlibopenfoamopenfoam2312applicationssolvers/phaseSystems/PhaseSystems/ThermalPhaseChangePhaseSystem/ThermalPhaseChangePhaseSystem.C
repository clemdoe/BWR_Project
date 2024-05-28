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

#include "ThermalPhaseChangePhaseSystem.H"
#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::iDmdt
(
    const phasePairKey& key
) const
{
    if (!iDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar dmdtSign(Pair<word>::compare(iDmdt_.find(key).key(), key));

    return dmdtSign**iDmdt_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::wDmdt
(
    const phasePairKey& key
) const
{
    if (!wDmdt_.found(key))
    {
        return phaseSystem::dmdt(key);
    }

    const scalar dmdtSign(Pair<word>::compare(wDmdt_.find(key).key(), key));

    return dmdtSign**wDmdt_[key];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
ThermalPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    volatile_(this->template lookupOrDefault<word>("volatile", "none")),
    saturationModel_
    (
        saturationModel::New(this->subDict("saturationModel"), mesh)
    ),
    phaseChange_(this->lookup("phaseChange"))
{

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        // Initially assume no mass transfer
        iDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("iDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        wDmdt_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("wDmdt", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        // Initially assume no mass transfer
        wMDotL_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("wMDotL", pair.name()),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
~ThermalPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::saturation() const
{
    return saturationModel_();
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    //Estos son purebas que fui haciendo
    //volScalarField idmdtDario (this->iDmdt(key));
//    volScalarField wdmdtDario (this->wDmdt(key));
//    volScalarField Basedmdt (BasePhaseSystem::dmdt(key));
//    //Foam::Info << "dmdt max: " << gMax(dmdtbis) << " dmdt min: " << gMin(dmdtbis) << Foam::endl;

//    Foam::Info<< " ----------------------------------------------- " << endl;
    //Foam::Info<< "Tasa transf. interfacial idmdtDario: " << idmdtDario.primitiveField() << endl;
//    Foam::Info<< "Tasa transf. pared: " << wdmdtDario.primitiveField() << endl;
//    Foam::Info<< "Base dmdt no se q es: " << Basedmdt.primitiveField() << endl;
//    Foam::Info<< " ----------------------------------------------- " << endl;

    return BasePhaseSystem::dmdt(key) + this->iDmdt(key) + this->wDmdt(key);
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(iDmdtTable, iDmdt_, iDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[iDmdtIter.key()];
        const volScalarField& iDmdt = *iDmdtIter();

        this->addField(pair.phase1(), "dmdt", iDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", - iDmdt, dmdts);
    }

    forAllConstIter(wDmdtTable, wDmdt_, wDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[wDmdtIter.key()];
        const volScalarField& wDmdt = *wDmdtIter();

        this->addField(pair.phase1(), "dmdt", wDmdt, dmdts);
        this->addField(pair.phase2(), "dmdt", - wDmdt, dmdts);
    }

    return dmdts;
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
// Termino de transferencia de energia que se agregar a la ecuacion de energia
// en el lado derecho de la ecuacion.
// Este creo q incorpora el termino fuente por transferncia de materia debio a la compresibilidad

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    //Foam::Info << " INGRESO A HEAT TRANSFER " << endl;

    // Add boundary term
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        if (this->wMDotL_.found(phasePairIter.key()))
        {
            const phasePair& pair(phasePairIter());

            if (pair.ordered())
            {
                continue;
            }

            const phaseModel& phase1 = pair.phase1();
            const phaseModel& phase2 = pair.phase2();

//            volScalarField phase1negpar (negPart(*this->wMDotL_[pair]));
//            volScalarField phase2negpar (posPart(*this->wMDotL_[pair]));
//            Foam::Info<< "Termino energia fase 1: " << phase1negpar.primitiveField() << endl;
//            Foam::Info<< "Termino energia fase 2: " << phase2negpar.primitiveField() << endl;

            *eqns[phase1.name()] += negPart(*this->wMDotL_[pair]);
            *eqns[phase2.name()] -= posPart(*this->wMDotL_[pair]);

            if
            (
                phase1.thermo().he().member() == "e"
             || phase2.thermo().he().member() == "e"
            )
            {
                const volScalarField dmdt
                (
                    this->iDmdt(pair) + this->wDmdt(pair)
                );

//                volScalarField Enphase1 (phase1.thermo().p()*dmdt/phase1.thermo().rho());
//                volScalarField Enphase2 (phase2.thermo().p()*dmdt/phase2.thermo().rho());
//                Foam::Info<< " ----------------------------------------------- " << endl;
//                Foam::Info<< "Termino energia fase 1: " << Enphase1.primitiveField() << endl;
//                Foam::Info<< "Termino energia fase 2: " << Enphase2.primitiveField() << endl;
//                Foam::Info<< " ----------------------------------------------- " << endl;

                if (phase1.thermo().he().member() == "e")
                {
                    *eqns[phase1.name()] +=
                        phase1.thermo().p()*dmdt/phase1.thermo().rho();
                }

                if (phase2.thermo().he().member() == "e")
                {
                    *eqns[phase2.name()] -=
                        phase2.thermo().p()*dmdt/phase2.thermo().rho();
                }
            }
        }
    }

    return eqnsPtr;
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::massTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::massTransfer() const
{
    autoPtr<phaseSystem::massTransferTable> eqnsPtr =
        BasePhaseSystem::massTransfer();

    phaseSystem::massTransferTable& eqns = eqnsPtr();

    Foam::Info << " INGRESO A MASS TRANSFER " << endl;

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        this->phasePairs_,
        phasePairIter
    )
    {
        const phasePair& pair(phasePairIter());

        if (pair.ordered())
        {
            continue;
        }

        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            if (Yi[i].member() != volatile_)
            {
                continue;
            }

            //Foam::Info << " INGRESO A MASS TRANSFER 2 " << endl;

            const word name
            (
                IOobject::groupName(volatile_, phase.name())
            );

            const word otherName
            (
                IOobject::groupName(volatile_, otherPhase.name())
            );

            // Note that the phase YiEqn does not contain a continuity error
            // term, so these additions represent the entire mass transfer

            const volScalarField dmdt(this->iDmdt(pair) + this->wDmdt(pair));

            Info << "<><<><><><><Entra en el calculo de las especies  y obtiene el dmdt <><><><><><" << endl;
            Info << "Coef de transf dmdt - max: " << max(dmdt.primitiveField()) << " min: " << min(dmdt.primitiveField())<< endl;


            *eqns[name] += dmdt;
            *eqns[otherName] -= dmdt;
        }
    }

    //Info << "Reportar en el massTransfer" << eqns << endl;

    return eqnsPtr;
}


template<class BasePhaseSystem>
void
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::correctInterfaceThermo()
{
    typedef compressible::alphatPhaseChangeWallFunctionFvPatchScalarField
        alphatPhaseChangeWallFunction;

    forAllConstIter
    (
        typename BasePhaseSystem::heatTransferModelTable,
        this->heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair
        (
            this->phasePairs_[heatTransferModelIter.key()]
        );

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

//        const volScalarField& T1(phase1.thermo().T());
//        const volScalarField& T2(phase2.thermo().T());

        // ------------------- N3---------------------------------- //
        //Esta es la modificacion N3

        volScalarField T1(phase1.thermo().T());
        volScalarField T2(phase2.thermo().T());

        const volScalarField& he1(phase1.thermo().he());
        const volScalarField& he2(phase2.thermo().he());

        const volScalarField& p(phase1.thermo().p());

        volScalarField& iDmdt(*this->iDmdt_[pair]);
        volScalarField& Tf(*this->Tf_[pair]);

        const volScalarField Tsat(saturationModel_->Tsat(phase1.thermo().p()));

//        Foam::Info<< " ----------------------------------------------- " << endl;
//        Info << "El nombre de la fase es: " << pair.name() << endl;

//        Info << "Temp. Fase 1: " <<  T1.primitiveField() << endl;
//        Info << "Temp. Fase 2: " << T2.primitiveField() << endl;
//        Info << "Entalpia h1: " << he1.primitiveField() << endl;
//        Info << "Entalpia h2: " << he2.primitiveField() << endl;
//        Info << "Presion p: " << p.primitiveField() << endl;
//        Info << "Transf. masa iDmdt: " << iDmdt.primitiveField() << endl;
//        Info << "Temp. Tf: " << Tf.primitiveField() << endl;
//        Info << "Temp. sat: " << Tsat.primitiveField() << endl;

//        Foam::Info<< " ----------------------------------------------- " << endl;

        volScalarField hf1
        (
            he1.member() == "e"
          ? phase1.thermo().he(p, Tsat) + p/phase1.rho()
          : phase1.thermo().he(p, Tsat)
        );
        volScalarField hf2
        (
            he2.member() == "e"
          ? phase2.thermo().he(p, Tsat) + p/phase2.rho()
          : phase2.thermo().he(p, Tsat)
        );

        volScalarField h1
        (
            he1.member() == "e"
          ? he1 + p/phase1.rho()
          : tmp<volScalarField>(he1)
        );

        volScalarField h2
        (
            he2.member() == "e"
          ? he2 + p/phase2.rho()
          : tmp<volScalarField>(he2)
        );

        volScalarField L
        (
            (neg0(iDmdt)*hf2 + pos(iDmdt)*h2)
          - (pos0(iDmdt)*hf1 + neg(iDmdt)*h1)
        );

//        volScalarField prueba1 (neg0(iDmdt)*hf2);
//        volScalarField prueba2 (pos(iDmdt)*h2);
//        volScalarField prueba3 (pos0(iDmdt)*hf1);
//        volScalarField prueba4 (neg(iDmdt)*h1);

        volScalarField iDmdtNew(iDmdt);

//        Foam::Info<< " ----------------------------------------------- " << endl;
//        Info << "Prueba 1: " << prueba1.primitiveField() << endl;
//        Info << "Prueba 2: " << prueba2.primitiveField() << endl;
//        Info << "Prueba 3: " << prueba3.primitiveField() << endl;
//        Info << "Prueba 4: " << prueba4.primitiveField() << endl;
//        Foam::Info<< " ----------------------------------------------- " << endl;


//        Foam::Info<< " ----------------------------------------------- " << endl;
//        Info << "El nombre de la fase es: " << pair.name() << endl;

//        Info << "Entalpia hf1: " <<  hf1.primitiveField() << endl;
//        Info << "Entalpia hf2: " << hf2.primitiveField() << endl;
//        Info << "Entalpia h1: " << h1.primitiveField() << endl;
//        Info << "Entalpia h2: " << h2.primitiveField() << endl;
//        Info << "Calor latente: " << L.primitiveField() << endl;
//        Info << "Transf. masa nuevo iDmdt: " << iDmdtNew.primitiveField() << endl;

//        Foam::Info<< " ----------------------------------------------- " << endl;


        // Ojo por esto lo agrego yo


        //Esto lo agregue en la Modificacion N1
//        volScalarField iDmdtLimit(iDmdt);
//        iDmdtLimit == (dimensionedScalar(iDmdt.dimensions(), 1e-6));



        // ------------------- N3---------------------------------- //
        //Modificacion N3 - Voy a limitar las temperaturas para que no
        //alla evaporacion en la fase gas

//        volScalarField Tsatur(Tf);

//        Tsatur == (dimensionedScalar(Tf.dimensions(), 530.65));

//        T1 = min(T1,Tsatur);
//        T2 = min(T2,Tsatur);


        if (phaseChange_)
        {
            volScalarField H1(heatTransferModelIter().first()->K(0));
            volScalarField H2(heatTransferModelIter().second()->K(0));

//            volScalarField H1sat (H1*(Tsat - T1));
//            volScalarField H2sat (H2*(Tsat - T2));

//            Foam::Info<< " ----------------------------------------------- " << endl;
//            Info << "Coef. de conveccion H1: " << H1.primitiveField() << endl;
//            Info << "Coef. de conveccion H2: " << H2.primitiveField() << endl;

//            Info << "Temp. Fase 1: " <<  T1.primitiveField() << endl;
//            Info << "Temp. Fase 2: " << T2.primitiveField() << endl;

//            Info << "H1sat" << H1sat.primitiveField() <<endl;
//            Info << "H2sat" << H2sat.primitiveField() <<endl;
//            Foam::Info<< " ----------------------------------------------- " << endl;

            //iDmdtNew = (H1*(Tsat - T1))/L;
            iDmdtNew = (H1*(Tsat - T1) + H2*(Tsat - T2))/L;


            // Esto lo agreugue yo. Es para limitar el valor a solo condensacion Modificacion N1
            //iDmdtNew = min(iDmdtNew,iDmdtLimit);

            //iDmdtNew = max(iDmdtNew,0.0);
        }
        else
        {
            iDmdtNew == dimensionedScalar(iDmdt.dimensions(), 0);
        }

        volScalarField H1(heatTransferModelIter().first()->K());
        volScalarField H2(heatTransferModelIter().second()->K());

        // Limit the H[12] to avoid /0
        H1.max(small);
        H2.max(small);

        //        Info << " --------------- Info extra ------------" << endl;
//        Info << "Coef de transf H1 - max: " << max(H1.primitiveField()) << " min: " << min(H1.primitiveField())<< endl;
//        Info << "Coef de transf H2 - max: " << max(H2.primitiveField()) << " min: " << min(H2.primitiveField())<< endl;
//        Info << "Coef de transf iDmdtNew - max: " << max(iDmdtNew.primitiveField()) << " min: " << min(iDmdtNew.primitiveField())<< endl;
//        Info << "Coef de transf iDmdt - max: " << max(iDmdt.primitiveField()) << " min: " << min(iDmdt.primitiveField())<< endl;
//        //Info << "Coef de transf dmdt - max: " << max(dmdt.primitiveField()) << " min: " << min(dmdt.primitiveField())<< endl;
//        Info << " ---------------------------------------" << endl;


        //Ojo porque aca estoy ponieod algo q no se si va a fucnionar. Le agrege el maximo
        //Esto lo agregue en la Modifiacion N1
//        volScalarField Tsatur(Tf);

//        Tsatur == (dimensionedScalar(Tf.dimensions(), 530.65));

        //Tf = min(Tsatur,(H1*T1 + H2*T2 + iDmdtNew*L)/(H1 + H2));

         //Tf = (H1*T1 + iDmdtNew*L)/(H1);

        //Esto lo agrego yo en la modificacion N1
//        iDmdt = min(iDmdtNew,iDmdtLimit);
//        Tf = min((H1*T1 + H2*T2 + iDmdtNew*L)/(H1 + H2),Tsatur);

        Tf = (H1*T1 + H2*T2 + iDmdtNew*L)/(H1 + H2);

        Info<< "Tf." << pair.name()
            << ": min = " << min(Tf.primitiveField())
            << ", mean = " << average(Tf.primitiveField())
            << ", max = " << max(Tf.primitiveField())
            << endl;

        scalar iDmdtRelax(this->mesh().fieldRelaxationFactor("iDmdt"));
        iDmdt = (1 - iDmdtRelax)*iDmdt + iDmdtRelax*iDmdtNew;

        // OJOOOOOOOOOOOOOOOO ESTO ES SOLO UNA PRUEBA!!!
        //Aca estoy haciendo cero la transferencia de masa por
        // por interface que no entiendo q es. ver la tesis
        // donde habla de esto

        //iDmdt = 1e-11*iDmdt;

        if (phaseChange_)
        {
            Info<< "iDmdt." << pair.name()
                << ": min = " << min(iDmdt.primitiveField())
                << ", mean = " << average(iDmdt.primitiveField())
                << ", max = " << max(iDmdt.primitiveField())
                << ", integral = " << fvc::domainIntegrate(iDmdt).value()
                << endl;
        }

        volScalarField& wDmdt(*this->wDmdt_[pair]);
        volScalarField& wMDotL(*this->wMDotL_[pair]);
        wDmdt = Zero;
        wMDotL = Zero;

        bool wallBoilingActive = false;

        forAllConstIter(phasePair, pair, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            if
            (
                phase.mesh().foundObject<volScalarField>
                (
                    "alphat." +  phase.name()
                )
            )
            {
                const volScalarField& alphat =
                    phase.mesh().lookupObject<volScalarField>
                    (
                        "alphat." +  phase.name()
                    );

                const fvPatchList& patches = this->mesh().boundary();
                forAll(patches, patchi)
                {
                    const fvPatch& currPatch = patches[patchi];

                    if
                    (
                        isA<alphatPhaseChangeWallFunction>
                        (
                            alphat.boundaryField()[patchi]
                        )
                    )
                    {
                        const alphatPhaseChangeWallFunction& PCpatch =
                            refCast<const alphatPhaseChangeWallFunction>
                            (
                                alphat.boundaryField()[patchi]
                            );

                        phasePairKey key(phase.name(), otherPhase.name());

                        if (PCpatch.activePhasePair(key))
                        {
                            wallBoilingActive = true;

                            const scalarField& patchDmdt =
                                PCpatch.dmdt(key);
                            const scalarField& patchMDotL =
                                PCpatch.mDotL(key);

                            const scalar sign
                            (
                                Pair<word>::compare(pair, key)
                            );

                            forAll(patchDmdt, facei)
                            {
                                const label faceCelli =
                                    currPatch.faceCells()[facei];
                                wDmdt[faceCelli] -= sign*patchDmdt[facei];
                                wMDotL[faceCelli] -= sign*patchMDotL[facei];
                            }
                        }
                    }
                }
            }
        }

        if (wallBoilingActive)
        {
            Info<< "wDmdt." << pair.name()
                << ": min = " << min(wDmdt.primitiveField())
                << ", mean = " << average(wDmdt.primitiveField())
                << ", max = " << max(wDmdt.primitiveField())
                << ", integral = " << fvc::domainIntegrate(wDmdt).value()
                << endl;
        }
    }
}


template<class BasePhaseSystem>
bool Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
