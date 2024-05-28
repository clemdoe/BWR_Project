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

#include "Grace.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Grace, 0);
    addToRunTimeSelectionTable(dragModel, Grace, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Grace::Grace
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
   dragModel(dict, pair, registerObject),
   residualRe_("residualRe", dimless, dict),
   muref1_("muref1", dimensionSet(1, -1, -1, 0, 0,0,0), dict),
   gref_("gref",dimensionSet(0,1,-2,0,0,0,0),dict)

{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Grace::~Grace()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::Grace::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(pair_.Eo());
    volScalarField Mo(pair_.Mo());


        //Valor de H de Grace que mesirve para luego evaluar J de Grace
        volScalarField Hgra
                (
                    (4/3) * Eo
                    *pow(Mo,-0.149)
                    *pow(((pair_.continuous().nu()*pair_.continuous().rho())/muref1_),-0.14)
                );

        //Valor de J de Grace
        volScalarField Jgra
                 (
                    neg(Hgra - 59.3)*0.94*pow(Hgra,0.757)
                    + pos(Hgra - 59.3)*3.42* pow(Hgra,0.441)
                 );


       //Velocidad terminal segun la ecuacion de Grace
       volScalarField UTer
                (
                pair_.continuous().nu()
                *pow(Mo,-0.149)
                *(Jgra-0.857)
                * (1.0/pair_.dispersed().d())
                );

    // --------------------------------------------------------------------------//
    //Esto solo lo uso cuando corro un caso sin gravedad
    //Esto lo hice cuando corri el caso de Marschall de flujo segregado que lo corria sin gravedad
    //Lo unico que cambia con el caso con gravedad es que el numero de Morton (Mo) en vez de llamarlo y
    //calcularlo con la gravedad del probllema lo calculo con una gravedad artificial y calculo el Moref

//    volScalarField Moref(mag(gref_)*pair_.continuous().nu()*pow3(pair_.continuous().nu()*pair_.continuous().rho()/pair_.sigma()));

//    //Valor de H de Grace que mesirve para luego evaluar J de Grace
//    volScalarField Hgra
//            (
//                (4/3) * Eo
//                *pow(Moref,-0.149)
//                *pow(((pair_.continuous().nu()*pair_.continuous().rho())/muref1_),-0.14)
//            );

//    //Info<< "BANANA" << endl;

//    //Valor de J de Grace
//    volScalarField Jgra
//             (
//                neg(Hgra - 59.3)*0.94*pow(Hgra,0.757)
//                + pos(Hgra - 59.3)*3.42* pow(Hgra,0.441)
//             );


//    volScalarField UTer
//             (
//             pair_.continuous().nu()
//             *pow(Moref,-0.149)
//             *(Jgra-0.857)
//             * (1.0/pair_.dispersed().d())
//             );

    // --------------------------------------------------------------------------//

   volScalarField Cdelip
           (
               (4/3)* gref_*pair_.dispersed().d()*
               (pair_.continuous().rho()-pair_.dispersed().rho())
               /(pair_.continuous().rho())
               *(1.0/sqr(UTer))
            );

   volScalarField Cdsphe
           (
              neg(Re - 0.01)*24/max(Re, residualRe_)
              + pos(Re - 0.01)*(24*(1+0.15*pow(Re,0.441)))/max(Re, residualRe_)
           );

   volScalarField Cdtot
           (
           max(Cdelip,Cdsphe)
           );

   // Ver el valor de la gravedad si lo pongo asi o de otra forma
    return
        Cdtot*max(Re, residualRe_);
}


// ************************************************************************* //
