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
    Foam::ddWallFluxMixedFvPatchScalarField

Description
    Abstract base class for scalar flux boundary conditions at solid surfaces,
    derived from mixedFvPatchScalarField.

    Enforces the particle flux at absorbing walls by setting the relation of 
    the face and the cell center value. The wall flux is imposed through the 
    valueFraction field of the mixed BC.

    This is an abstract class. It cannot be used directly in a case file
    as the 'type'. Derived classes must implement the pure virtual function.

Usage
    This class is not instantiated directly. Derived classes are instantiated.
    
SourceFiles
    ddWallFluxMixedFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef ddWallFluxMixedFvPatchScalarField_H
#define ddWallFluxMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "SoPLASMAConstants.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "driftDiffusion.H"
#include "plasmaWallBC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class ddWallFluxMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ddWallFluxMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public plasmaWallBC
{
protected:

    // Protected Data

        //- The name of the temperature field (e.g., "T_e"). 
        word TName_;

        //- Constant Temperature Value
        dimensionedScalar TValue_;

    // Protected Member Functions

        //- Calculate thermal velocity (account for the 1/4)
        //- For a single constant temperature (dimensionedScalar)
        dimensionedScalar calcThermalVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T
        ) const;

        //- Calculate thermal velocity (account for the 1/4)
        //- For a field temperature (scalarField)
        tmp<scalarField> calcThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
        ) const;

        //- Calculate the wall absorption velocity
        virtual tmp<scalarField> calcAbsorptionVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const = 0;

        //- Calculate the effective wall velocity
        virtual tmp<scalarField> calcEffectiveWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal
        ) const = 0;

public:

    //- Runtime type information
    TypeName("ddWallFluxMixed");

    // Constructors

        //- Construct from patch and internal field
        ddWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ddWallFluxMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        ddWallFluxMixedFvPatchScalarField
        (
            const ddWallFluxMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ddWallFluxMixedFvPatchScalarField
        (
            const ddWallFluxMixedFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ddWallFluxMixedFvPatchScalarField
        (
            const ddWallFluxMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const = 0;

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const = 0;

        // Member Functions

            //- Update the coefficients
            virtual void updateCoeffs();

            //- Returns the enableSurfaceCharging flag
            virtual bool enableSurfaceCharging() const override
            {
                return false; 
            }

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
