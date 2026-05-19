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
    Foam::neutralDDWallFluxMixedFvPatchScalarField

Description
    Concrete implementation of ddWallFluxMixedFvPatchScalarField for 
    neutral species.

    This boundary condition calculates the neutral flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            neutralDDWallFluxMixed;

        // Neutral temperature [K]
        T               Tn; 
    }
    \endverbatim
    
SourceFiles
    neutralDDWallFluxMixedFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef neutralDDWallFluxMixedFvPatchScalarField_H
#define neutralDDWallFluxMixedFvPatchScalarField_H

#include "ddWallFluxMixedFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class neutralDDWallFluxMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class neutralDDWallFluxMixedFvPatchScalarField
:
    public ddWallFluxMixedFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the neutral-specific wall absorption velocity
        virtual tmp<scalarField> calcAbsorptionVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

        //- Calculate the neutral-specific wall effective velocity
        virtual tmp<scalarField> calcEffectiveWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

public:

    //- Runtime type information
    TypeName("neutralDDWallFluxMixed");

    // Constructors

        //- Construct from patch and internal field
        neutralDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neutralDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        neutralDDWallFluxMixedFvPatchScalarField
        (
            const neutralDDWallFluxMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        neutralDDWallFluxMixedFvPatchScalarField
        (
            const neutralDDWallFluxMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        neutralDDWallFluxMixedFvPatchScalarField
        (
            const neutralDDWallFluxMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new neutralDDWallFluxMixedFvPatchScalarField(*this)
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
                new neutralDDWallFluxMixedFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
