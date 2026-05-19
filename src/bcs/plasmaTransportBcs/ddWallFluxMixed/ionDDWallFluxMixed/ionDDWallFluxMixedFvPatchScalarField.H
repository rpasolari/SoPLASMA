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
    Foam::ionDDWallFluxMixedFvPatchScalarField

Description
    Concrete implementation of ddWallFluxMixedFvPatchScalarField for 
    ion species.

    This boundary condition calculates the ion flux to a solid surface 
    based on thermal motion. The drift flux is also added, if it pushes ions
    towards the boundary.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            ionDDWallFluxMixed;

        // Ion temperature [K]
        T               Tion; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging false;
    }
    \endverbatim
    
SourceFiles
    ionDDWallFluxMixedFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef ionDDWallFluxMixedFvPatchScalarField_H
#define ionDDWallFluxMixedFvPatchScalarField_H

#include "ddWallFluxMixedFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class ionDDWallFluxMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ionDDWallFluxMixedFvPatchScalarField
:
    public ddWallFluxMixedFvPatchScalarField
{
protected:

    // Protected data

        //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

    // Protected Member Functions

        //- Calculate the ion-specific wall absorption velocity
        virtual tmp<scalarField> calcAbsorptionVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

        //- Calculate the ion-specific wall effective velocity
        virtual tmp<scalarField> calcEffectiveWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const override;

public:

    //- Runtime type information
    TypeName("ionDDWallFluxMixed");

    // Constructors

        //- Construct from patch and internal field
        ionDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ionDDWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        ionDDWallFluxMixedFvPatchScalarField
        (
            const ionDDWallFluxMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ionDDWallFluxMixedFvPatchScalarField
        (
            const ionDDWallFluxMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ionDDWallFluxMixedFvPatchScalarField
        (
            const ionDDWallFluxMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new ionDDWallFluxMixedFvPatchScalarField(*this)
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
                new ionDDWallFluxMixedFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Return the surface charging switch value
            inline bool enableSurfaceCharging() const
            {
                return enableSurfaceCharging_;
            }
            
            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
