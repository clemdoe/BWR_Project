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


Class
    Foam::dragModels::GraceCorrelated

Description
    Correlation of Grace et al.

    Reference:
    \verbatim
        Clift, R., Grace, J. R., & Weber, M. E. (2005). Bubbles, drops, and particles. Courier Corporation.
    \endverbatim

    // Las expresion utilizadas estan en el Capitulo 7 - Elipsoidal Fluid Particles

SourceFiles
    GraceCorrelated.C



\*---------------------------------------------------------------------------*/


#ifndef Grace_H
#define Grace_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class GraceCorrelated Declaration
\*---------------------------------------------------------------------------*/

class Grace
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;
        //- Viscosidad de referencia para agua a 25 grados. Libro de Grace pagina 176.
        //- Se lee desde diccionario y vale 0.0009 kg / ms
        const dimensionedScalar muref1_;
        //- Aceleracion de la gravedad. Aca la defino a pata no mas. Las unidades estan
        //- definidas en Grace.C
        const dimensionedScalar gref_;

public:

    //- Runtime type information
    TypeName("Grace");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Grace
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Grace();

    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
