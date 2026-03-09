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
    Foam::tabulatedMobility1D

Description
    1D tabulated mobility model for plasma species.

    This model calculates the mobility mu as a function of a user-defined
    variable, by interpolating from a user-provided data file. The model
    supports both native OpenFOAM list format and CSV format.

    Usage - Native OpenFOAM Format:
    \verbatim
        mobilityCoeffs
        {
            lookupVariable  E;           // Mandatory: Name of field for lookup
            file            "path";      // Mandatory: Data file path
            readerType      openFoam;    // Optional:  Default is openFoam
            outOfBounds     clamp;       // Optional:  Default is warn
        }
    \endverbatim

    Usage - CSV Format:
    \verbatim
        mobilityCoeffs
        {
            lookupVariable      E;       // Mandatory
            readerType          csv;     // Mandatory
            file                "path";  // Mandatory
            hasHeaderLine       true;    // Mandatory: true to skip first line
            refColumn           0;       // Mandatory: lookup column index
            componentColumns    (1);     // Mandatory: value column index list
            separator           ",";     // Optional:  Default is ","
            outOfBounds         clamp;   // Optional:  Default is warn
        }
    \endverbatim

    The 'outOfBounds' entry handles values outside the table range:
    - clamp:        Use the value of the nearest endpoint (Recommended).
    - error:        Stop the simulation with a Fatal Error.
    - warn:         Clamp the value but print a warning to the log.
    - repeat:       Treat the table as a periodic cycle.

    The native data file should contain a list of (fieldValues mu_value) pairs:
    \verbatim
        (
            (0.0      1.2e-4)
            (100.0    1.1e-4)
            ...
        )
    \endverbatim

    It is selected at runtime via:
    \verbatim
        mobilityModel   tabulatedMobility1D;
    \endverbatim

SourceFiles
    tabulatedMobility1D.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tabulatedMobility1D_H
#define Foam_tabulatedMobility1D_H

#include "interpolationTable.H"

#include "plasmaMobilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedMobility1D Declaration
\*---------------------------------------------------------------------------*/

class tabulatedMobility1D
:
    public plasmaMobilityModel
{
private:

    //- The 1D interpolation table (scalar -> scalar)
    interpolationTable<scalar> table_;

    //- Variable to use for lookup: "E", "EN", "Te" etc
    word lookupVariable_;

public:

    //- Runtime type information
    TypeName("tabulated1D");

    // Constructors

        //- Construct from components
        tabulatedMobility1D
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
