/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
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
    Foam::plasmaMobilityModel

Description
    Physics manager for plasma species mobility modeling.

Usage
    The model is selected via the 'type' entry in the mobility dictionary
    within the species settings.

    \verbatim
    mobility
    {
        // 1. Constant Mobility
        type            constant;
        value           1.5e-4;      // Mandatory

        // 2. Tabulated 1D Mobility
        // mu = f(x) interpolated from a data file.
        type            tabulated1D;
        
        // Native OpenFOAM Format usage:
        lookupVariable  E;           // Mandatory
        file            "path";      // Mandatory
        readerType      openFoam;    // Optional:  Default is openFoam
        outOfBounds     clamp;       // Optional:  Default is warn

        // CSV Format usage:
        readerType      csv;         // Mandatory
        hasHeaderLine   true;        // Mandatory: skip first line
        refColumn       0;           // Mandatory: lookup column index
        componentColumns (1);        // Mandatory: value column index list
        separator       ",";         // Optional:  Default is ","
        outOfBounds     clamp;       // Optional:  Default is warn

        // 3. Tabulated 2D Mobility
        // mu = f(x, y) via bilinear interpolation from a user data file.
        // File contains a nested list: ( X ( (Y1 mu1) (Y2 mu2) ) )
        type            tabulated2D;
        lookupVariables (E Te);      // Mandatory: List of two field names
        file            "path";      // Mandatory: Data file path
        outOfBounds     clamp;       // Optional: (clamp|error|warn)

        // 4. Power Law Mobility (Analytic)
        // mu = amplitude * |lookupVariable|^exponent
        type            powerLaw;
        amplitude       1.5e-4;      // Mandatory: muCoeff
        exponent        -0.5;        // Mandatory: fieldExp
        lookupVariable  E;           // Optional:  Default is E

        // 5. Universal Function1 Property
        // Supports any Function1 type: coded, polynomial, sine, etc.
        type            function1;
        lookupVariable  E;           // Mandatory
        function        coded;       // Function1 type (coded|polynomial|...)
        
        // Example for 'coded' type:
        name            myCustomMu;
        code            #{ return 1.5e-4 * Foam::exp(-x/100.0); #};

        // Example for 'polynomial' type:
        // function     polynomial ((1 0) (2 1)); // mu = 1 + 2*E
    }
    \endverbatim

    Notes on Out-Of-Bounds handling:
    - clamp:  Use value of the nearest endpoint (Recommended).
    - error:  Stop the simulation with a Fatal Error.
    - warn:   Clamp the value but print a warning to the log (Default).
    - repeat: Treat the table as a periodic cycle (1D only).

SourceFiles
    plasmaMobilityModel.C

\*---------------------------------------------------------------------------*/


#ifndef plasmaMobilityModel_H
#define plasmaMobilityModel_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "plasmaSpecies.H"
#include "plasmaPropertyEvaluator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class plasmaMobilityModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaMobilityModel
{
protected:

    //- The specific math evaluator (e.g., tabulated1D, powerLaw, etc.)
    autoPtr<plasmaPropertyEvaluator> evaluator_;

    //- Model type (e.g. "constant", "tabulated1D", etc.)
    const word modelType_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    const plasmaSpecies& species_;

    //- Index of the species this model applies to
    const label specieIndex_;

    //- Mobility model dictionary
    const dictionary dict_;

    //- Flag for uniform or non-uniform field property
    bool isUniform_;

    //- Flag for constant or dynamic field property
    bool isConstant_;

    //- Mutable flag to check if a constant field has already been calculated
    mutable bool upToDate_;

    //- Storage for the uniform value
    mutable dimensionedScalar muValue_;

    //- Storage for the volScalarField
    mutable autoPtr<volScalarField> muFieldPtr_;

    //- Private update logic
    void update() const;
    
public:

    //- Runtime type information
    TypeName("plasmaMobilityModel");

    // Constructors

        //- Construct base model
        plasmaMobilityModel
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex
        );

    //- Selector
    static autoPtr<plasmaMobilityModel> New
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    );

    //- Destructor
    virtual ~plasmaMobilityModel() = default;

    // Member Functions

        //- Returns true if the mobility is spatially uniform
        bool isUniform() const { return isUniform_; }

        //- Returns the actual scalar value (if spatially uniform)
        const dimensionedScalar& muValue() const { return muValue_; }

        //- Return mobility as a temporary field
        tmp<volScalarField> mu() const;

        //- Correct mobility
        void correct(volScalarField& mu) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
