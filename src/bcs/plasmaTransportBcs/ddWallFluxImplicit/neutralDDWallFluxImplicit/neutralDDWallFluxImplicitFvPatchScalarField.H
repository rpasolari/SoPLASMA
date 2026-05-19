/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::neutralDDWallFluxImplicitFvPatchScalarField

Description
    Concrete implementation of ddWallFluxImplicitFvPatchScalarField for 
    neutral species.

    This boundary condition calculates the neutral flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            neutralDDWallFluxImplicit;

        // Neutral temperature [K]
        T               Tn; 
    }
    \endverbatim
    
SourceFiles
    neutralDDWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef neutralDDWallFluxImplicitFvPatchScalarField_H
#define neutralDDWallFluxImplicitFvPatchScalarField_H

#include "ddWallFluxImplicitFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class neutralDDWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class neutralDDWallFluxImplicitFvPatchScalarField
:
    public ddWallFluxImplicitFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Return u_th/4 for neutrals
        virtual tmp<scalarField> calcWallThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
        ) const override;

public:

    //- Runtime type information
    TypeName("neutralDDWallFluxImplicit");

    // Constructors

        //- Construct from patch and internal field
        neutralDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neutralDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        neutralDDWallFluxImplicitFvPatchScalarField
        (
            const neutralDDWallFluxImplicitFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        neutralDDWallFluxImplicitFvPatchScalarField
        (
            const neutralDDWallFluxImplicitFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        neutralDDWallFluxImplicitFvPatchScalarField
        (
            const neutralDDWallFluxImplicitFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new neutralDDWallFluxImplicitFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new neutralDDWallFluxImplicitFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- No drift for neutrals: always return zero
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const override;

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
