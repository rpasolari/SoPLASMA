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
    Foam::ddWallFluxImplicitFvPatchScalarField

Description
    Abstract base class for wall flux boundary conditions in drift-diffusion
    plasma transport that acts directly on the linear system matrix.

    Rather than setting refValue and refGradient explicitly, this class
    modifies the diagonal and source coefficients of the fvScalarMatrix
    assembled for the species transport equation. The wall flux, composed
    of the thermal flux and the drift contribution, is incorporated as an
    implicit source term, coupling the boundary cell value directly into
    the matrix.

    This is an abstract class. It cannot be used directly in a case file
    as the 'type'. Derived classes must implement the pure virtual function.

Usage
    This class is not instantiated directly. Derived classes are instantiated.
    
SourceFiles
    ddWallFluxImplicitFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_ddWallFluxImplicitFvPatchScalarField_H
#define Foam_ddWallFluxImplicitFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fvPatchField.H"
#include "fvPatchFieldMapper.H"

#include "SoPLASMAConstants.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "driftDiffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class ddWallFluxImplicitFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ddWallFluxImplicitFvPatchScalarField
:
    public fvPatchScalarField
{
protected:

    // Protected Data

        //- The name of the temperature field (e.g., "T_e"). 
        word TName_;

        //- Constant Temperature Value
        dimensionedScalar TValue_;

    // Protected Member Functions

        //- Calculate thermal velocity u_th/4 for a constant temperature
        dimensionedScalar calcThermalVelocity
        (
            const dimensionedScalar& m,
            const dimensionedScalar& T
        ) const;

        //- Calculate thermal velocity u_th/4 for a temperature field
        tmp<scalarField> calcThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
        ) const;

        //- Pure virtual: return u_th/4 for this species on this patch
        //  Derived classes implement this based on species mass and temperature
        virtual tmp<scalarField> calcWallThermalVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T
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

            //- Update the face value (n_face = n_cell, consistent with matrix)
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Update the coefficients
            virtual void updateCoeffs();

            //- Coefficient for fvm::div diagonal
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- RHS coefficient for fvm::div
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Coefficient for fvm::laplacian diagonal
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- RHS coefficient for fvm::laplacian
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

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
