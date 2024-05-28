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

Class
    Foam::temperatureCoupledBase_Modif2

Description

    Es igual al original pero lo modifique para poder utilizar la condicion de borde
    turbulentTemperatureCoupledBaffleMixed en el problema de dos fluidos

    Common functions used in temperature coupled boundaries.

    The thermal conductivity \c kappa may be obtained by the following methods:
      - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
        defined by 'kappa'
      - 'fluidThermo' : use fluidThermo and default
        compressible::turbulenceModel to calculate kappa
      - 'solidThermo' : use solidThermo kappa()
      - 'directionalSolidThermo': uses look up for volSymmTensorField for
        transformed kappa vector. Field name definable in 'alphaAni',
        named 'Anialpha' in solid solver by default

    \par Keywords provided by this class:
      \table
        Property     | Description                | Required    | Default value
        kappaMethod  | Thermal conductivity method        | yes |
        kappa        | Name of thermal conductivity field | no  | none
        alphaAni     | Name of the non-isotropic alpha    | no  | Anialpha
      \endtable

Usage
    \verbatim
    nonIsotropicWall
    {
        ...
        kappaMethod     directionalSolidThermo;
        kappa           none;
        alphaAni        Anialpha;
        ...
    }
    \endverbatim

SourceFiles
    temperatureCoupledBase_Modif2.C

\*--------------------------_Modif2-------------------------------------------------*/

#ifndef temperatureCoupledBase_Modif2_H
#define temperatureCoupledBase_Modif2_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"


// ------ Esto lo incluyo para leer las propiedades de los dos fluidos ---- //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class temperatureCoupledBase_Modif2 Declaration
\*---------------------------------------------------------------------------*/

class temperatureCoupledBase_Modif2
{
public:

    // Public enumerations

        //- Type of supplied Kappa
        enum KMethodType
        {
            mtFluidThermo,
            mtSolidThermo,
            mtDirectionalSolidThermo,
            mtLookup,
            TwoFluid
        };


protected:

    // Protected data

        static const NamedEnum<KMethodType, 5> KMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- How to get K
        const KMethodType method_;

        //- Name of thermal conductivity field (if looked up from database)
        const word kappaName_;

        //- Name of the non-Isotropic alpha (default: Anialpha)
        const word alphaAniName_;


public:

    // Constructors

        //- Construct from patch and K name
        temperatureCoupledBase_Modif2
        (
            const fvPatch& patch,
            const word& calculationMethod,
            const word& kappaName,
            const word& alphaAniName
        );

        //- Construct from patch and dictionary
        temperatureCoupledBase_Modif2
        (
            const fvPatch& patch,
            const dictionary& dict
        );

         //- Construct from patch and temperatureCoupledBase_Modif2
        temperatureCoupledBase_Modif2
        (
            const fvPatch& patch,
            const temperatureCoupledBase_Modif2& base
        );


    // Member Functions

        //- Method to obtain K
        word KMethod() const
        {
            return KMethodTypeNames_[method_];
        }

        //- Name of thermal conductivity field
        const word& kappaName() const
        {
            return kappaName_;
        }

        //- Given patch temperature calculate corresponding K field
        tmp<scalarField> kappa(const scalarField& Tp) const;


        //- Creo un escalar cualquiera
        //tmp<scalarField> tw(new scalarField(patch_.size(), 1.1));

        //tmp<scalarField> tw(new scalarField(patch_.size(), 1.1));

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
