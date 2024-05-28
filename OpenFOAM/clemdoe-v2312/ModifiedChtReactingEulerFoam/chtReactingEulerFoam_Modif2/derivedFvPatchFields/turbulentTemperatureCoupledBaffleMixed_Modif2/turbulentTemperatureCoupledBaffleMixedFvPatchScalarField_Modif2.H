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
    Foam::compressible::
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles. Optional thin thermal layer resistances can be
    specified through thicknessLayers and kappaLayers entries.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = zero gradient
      - refValue = neighbour value
      - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    where KDelta is heat-transfer coefficient K * deltaCoeffs

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase_Modif2.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivities per layer [W/m/K] | no |
        kappaMethod  | inherited from temperatureCoupledBase_Modif1 | inherited |
        kappa        | inherited from temperatureCoupledBase_Modif1 | inherited |
        TwoFluidAcoupled | Bandera para indicar si en ese patch esta conectado con una region de dos fluidos | no

    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        Tnbr            T;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    // ---------------------------------------------- //

    TwoFluidAcoupled: ES un bool que si esto conectado una region de dos-fluidos con un
    solido o con una region monofasico pongo on y sino en caso de conectar un fluido monofasico con
    un solido va off

    // ---------------------------------------------- //

See also
    Foam::temperatureCoupledBase_Modif2

SourceFiles
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif1.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2_H
#define turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase_Modif2.H"
#include "scalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class turbulentTemperatureCoupledBaffleMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase_Modif2
{
    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;

        //- Bandera para indicar si en ese patch esta conectado con una
        //- region de dos fluidos
        bool TwoFluidAcoupled_;
        //const word TwoFluidAcoupled_;

        //- Name del tipo de region
        //- Esto se setea asi:
        //- Si estas en el T.water -> TypeRegion TwoFluid
        //- Si estas en el T acoplado a dos fluidos -> TypeRegion solido
        //- Si estas en el T acoplado pero a fluidos monofasicos -> TypeRegion fluido
        const word TypeRegion_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledBaffleMixed_Modif2");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
        (
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
        (
            const turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_Modif2
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
