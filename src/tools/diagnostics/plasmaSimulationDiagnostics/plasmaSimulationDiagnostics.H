/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
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
    Foam::plasmaSimulationDiagnostics

Description
    Class for calculating and reporting runtime diagnostics in plasma
    simulations.

    This class reads the 'plasmaSimulationDiagnostics' dictionary to determine
    which fields (e.g., species number densities, electric field, surface
    charge) should have their minimum and maximum values reported to the console
    at each time step.

SourceFiles
    plasmaSimulationDiagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaSimulationDiagnostics_H
#define plasmaSimulationDiagnostics_H

#include "fvCFD.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class plasmaTransport;
/*---------------------------------------------------------------------------*\
                 Class plasmaSimulationDiagnostics Declaration
\*---------------------------------------------------------------------------*/

class plasmaSimulationDiagnostics
{
private:
// Private data

        //- Reference to the Time object
        const Time& runTime_;

        //- Reference to the plasma transport model
        const plasmaTransport& transport_;

        //- Dictionary of control parameters
        dictionary dict_;

    // Diagnostic switches

        //- Switch to enable/disable printing species min/max values
        Switch printSpecies_;

        //- Switch to enable/disable printing electromagnetic fields min/max
        Switch printElectromagnetics_;

    // Private Member Functions

        //- Read controls from dictionary
        void read();

public:

    // Constructors

        //- Construct from Time and plasmaTransport
        plasmaSimulationDiagnostics
        (
            const Time& runTime,
            const plasmaTransport& transport
        );

    //- Destructor
    virtual ~plasmaSimulationDiagnostics() {}

    // Member Functions

        //- Report the requested diagnostics
        void report();

    // Accessors (Getters)

        //- Is species reporting enabled?
        inline bool printSpecies() const
        {
            return printSpecies_;
        }

        //- Is electromagnetics reporting enabled?
        inline bool printElectromagnetics() const
        {
            return printElectromagnetics_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
