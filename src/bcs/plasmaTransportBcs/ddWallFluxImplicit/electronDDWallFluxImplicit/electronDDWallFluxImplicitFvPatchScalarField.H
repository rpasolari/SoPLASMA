/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

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
    Foam::electronDDWallFluxImplicitFvPatchScalarField

Description
    Concrete implementation of ddWallFluxImplicitFvPatchScalarField for 
    electron species.

    This boundary condition calculates the electron flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            electronDDWallFluxImplicit;

        // Electron temperature [K]
        T               Te; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging_ false;
    }
    \endverbatim
    
SourceFiles
    electronDDWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef electronDDWallFluxImplicitFvPatchScalarField_H
#define electronDDWallFluxImplicitFvPatchScalarField_H

#include "ddWallFluxImplicitFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class electronDDWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class electronDDWallFluxImplicitFvPatchScalarField
:
    public ddWallFluxImplicitFvPatchScalarField
{
protected:

    // Protected data

       //- Switch to enable surface charging
        Switch enableSurfaceCharging_;

        //- Switch to include drift velocity in wall flux calculation
        Switch includeDriftFlux_;

        //- Switch to enable secondary electron emission
        Switch enableSEE_;

        //- Default SEE coefficient
        scalar defaultSEEC_;

        //- Species-specific SEE coefficients
        dictionary speciesSEEC_;

        //- List to keep the SEE coefficients for each ion species
        List<scalar> seec_;

        //- Switch to use refValue for the SEE instead of refGrad
        Switch setRefValue_;

        //- Switch to read the SEE entries only once
        Switch mapped_;

    // Protected Member Functions

        //- Calculate the electron-specific wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal,
            const scalar Z
        ) const override;

public:

    //- Runtime type information
    TypeName("electronDDWallFluxImplicit");

    // Constructors

        //- Construct from patch and internal field
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        electronDDWallFluxImplicitFvPatchScalarField
        (
            const electronDDWallFluxImplicitFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new electronDDWallFluxImplicitFvPatchScalarField(*this)
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
                new electronDDWallFluxImplicitFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            virtual void updateCoeffs();

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
