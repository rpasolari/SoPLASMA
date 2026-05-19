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
    Foam::ionDDWallFluxImplicitFvPatchScalarField

Description
    Concrete implementation of ddWallFluxImplicitFvPatchScalarField for 
    ion species.

    This boundary condition calculates the ion flux to a solid surface 
    based on thermal motion. The drift flux is also added, if it pushes ions
    towards the boundary

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            ionDDWallFluxImplicit;

        // Ion temperature [K]
        T               Tion; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging_ false;
    }
    \endverbatim
    
SourceFiles
    ionDDWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef ionDDWallFluxImplicitFvPatchScalarField_H
#define ionDDWallFluxImplicitFvPatchScalarField_H

#include "ddWallFluxImplicitFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class ionDDWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ionDDWallFluxImplicitFvPatchScalarField
:
    public ddWallFluxImplicitFvPatchScalarField
{
protected:

    // Protected data

        //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

    // Protected Member Functions

        //- Return u_th/4 for ions
        //  Simply delegates to the base calcThermalVelocity
        virtual tmp<scalarField> calcWallThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
        ) const override;

public:

    //- Runtime type information
    TypeName("ionDDWallFluxImplicit");

    // Constructors

        //- Construct from patch and internal field
        ionDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ionDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        ionDDWallFluxImplicitFvPatchScalarField
        (
            const ionDDWallFluxImplicitFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ionDDWallFluxImplicitFvPatchScalarField
        (
            const ionDDWallFluxImplicitFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ionDDWallFluxImplicitFvPatchScalarField
        (
            const ionDDWallFluxImplicitFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new ionDDWallFluxImplicitFvPatchScalarField(*this)
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
                new ionDDWallFluxImplicitFvPatchScalarField(*this, iF)
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
