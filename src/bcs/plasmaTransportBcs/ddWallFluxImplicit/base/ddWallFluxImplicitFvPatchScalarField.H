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
    Foam::ddWallFluxImplicitFvPatchScalarField

Description
    Abstract base class for scalar flux boundary conditions at solid surfaces,
    derived from mixedFvPatchScalarField.

    This class serves as a foundation for specific particle or charge flux 
    conditions (e.g., ion or electron fluxes) where particles impact a wall.

    This is an abstract class. It cannot be used directly in a case file
    as the 'type'. Derived classes must implement the pure virtual function.

Usage
    This class is not instantiated directly. Derived classes are instantiated.
    
SourceFiles
    ddWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_ddWallFluxImplicitFvPatchScalarField_H
#define Foam_ddWallFluxImplicitFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of plasmaTransport and plasmaSpecies classes
class plasmaTransport;
class plasmaSpecies;

/*---------------------------------------------------------------------------*\
         Class ddWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ddWallFluxImplicitFvPatchScalarField
:
    public mixedFvPatchScalarField
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

        //- Calculate the effective wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal,
            const scalar Z
        ) const = 0;

public:

    //- Runtime type information
    TypeName("ddWallFluxImplicit");

    // Constructors

        //- Construct from patch and internal field
        ddWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ddWallFluxImplicitFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        ddWallFluxImplicitFvPatchScalarField
        (
            const ddWallFluxImplicitFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ddWallFluxImplicitFvPatchScalarField
        (
            const ddWallFluxImplicitFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        ddWallFluxImplicitFvPatchScalarField
        (
            const ddWallFluxImplicitFvPatchScalarField&,
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
            virtual bool enableSurfaceCharging() const
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
