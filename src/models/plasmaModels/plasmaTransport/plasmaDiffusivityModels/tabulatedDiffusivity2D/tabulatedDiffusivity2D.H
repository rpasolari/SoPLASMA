/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::tabulatedDiffusivity2D

Description
    2D tabulated diffusivity model for plasma species.

    This model calculates the diffusivity D as a function of two user-defined
    variables (e.g., reduced electric field E/N and gas temperature T), 
    by performing a bilinear interpolation from a user-provided data file.

    Usage:
    \verbatim
        diffusivityCoeffs
        {
            lookupVariables  (E T);       // Mandatory: List of two field names
            file             "path";      // Mandatory: Data file path
            outOfBounds      clamp;       // Optional:  (clamp|error|warn)
        }
    \endverbatim

    The 'outOfBounds' entry handles values outside the table range:
    - clamp:  Use the value of the nearest endpoint (Recommended).
    - error:  Stop the simulation with a Fatal Error.
    - warn:   Clamp the value but print a warning to the log (Default).

    Native data file format:
    The file should contain a nested list where the first value is the X-axis 
    coordinate (e.g., E), followed by a list of (Y D) pairs:
    \verbatim
        (
            0.0 ( 
                (300.0  1.2e-4) 
                (500.0  1.1e-4) 
            )
            100.0 (
                (300.0  1.0e-4)
                (500.0  0.9e-4)
            )
        )
    \endverbatim

    It is selected at runtime via:
    \verbatim
        diffusivityModel   tabulatedDiffusivity2D;
    \endverbatim

SourceFiles
    tabulatedDiffusivity2D.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tabulatedDiffusivity2D_H
#define Foam_tabulatedDiffusivity2D_H

#include "interpolation2DTable.H"

#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedDiffusivity2D Declaration
\*---------------------------------------------------------------------------*/

class tabulatedDiffusivity2D
:
    public plasmaDiffusivityModel
{
private:

    //- The 2D interpolation table (scalar, scalar -> scalar)
    interpolation2DTable<scalar> table_;

    //- Variable to use for lookup: "E", "EN", "Te" etc
    wordList lookupVariables_;

public:

    //- Runtime type information
    TypeName("tabulated2D");

    // Constructors

        //- Construct from components
        tabulatedDiffusivity2D
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    // Member Functions
        
        //- Correct model
        virtual void correct(volScalarField& mu) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
