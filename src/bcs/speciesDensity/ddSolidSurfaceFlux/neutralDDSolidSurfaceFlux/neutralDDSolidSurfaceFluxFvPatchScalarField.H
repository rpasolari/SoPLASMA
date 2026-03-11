/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
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
    Foam::neutralDDSolidSurfaceFluxFvPatchScalarField

Description
    Concrete implementation of ddSolidSurfaceFluxFvPatchScalarField for 
    neutral species.

    This boundary condition calculates the neutral flux to a solid surface 
    based on thermal motion.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
        type            neutralDDSolidSurfaceFlux;

        // Electron temperature [K]
        T               Tn; 

        // Optional: Enable surface charge accumulation
        enableSurfaceCharging_ false;
    }
    \endverbatim
    
SourceFiles
    neutralDDSolidSurfaceFluxFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef neutralDDSolidSurfaceFluxFvPatchScalarField_H
#define neutralDDSolidSurfaceFluxFvPatchScalarField_H

#include "ddSolidSurfaceFluxFvPatchScalarField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class neutralDDSolidSurfaceFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class neutralDDSolidSurfaceFluxFvPatchScalarField
:
    public ddSolidSurfaceFluxFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the neutral-specific wall absorption velocity
        virtual tmp<scalarField> calcWallVelocity
        (
            const dimensionedScalar& m,
            const scalarField& T,
            const scalarField& uDriftNormal,
            const scalar Z
        ) const override;

public:

    //- Runtime type information
    TypeName("neutralDDSolidSurfaceFlux");

    // Constructors

        //- Construct from patch and internal field
        neutralDDSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neutralDDSolidSurfaceFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BC onto a new patch
        neutralDDSolidSurfaceFluxFvPatchScalarField
        (
            const neutralDDSolidSurfaceFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        neutralDDSolidSurfaceFluxFvPatchScalarField
        (
            const neutralDDSolidSurfaceFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        neutralDDSolidSurfaceFluxFvPatchScalarField
        (
            const neutralDDSolidSurfaceFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new neutralDDSolidSurfaceFluxFvPatchScalarField(*this)
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
                new neutralDDSolidSurfaceFluxFvPatchScalarField(*this, iF)
            );
        }

        // Member Functions

            //- Write dictionary entries
            virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
