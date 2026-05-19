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
    Foam::plasmaTimeControl

Description
    Class for automatic time-step control in plasma simulations.

    This class serves as a specialized replacement for the standard OpenFOAM
    'setDeltaT.H' functionality. It calculates the maximum permissible time
    step based on plasma-specific stability criteria.

    The class reads the 'plasmaTimeControl' dictionary (found within
    system/plasmaSimulationControl) to determine which criteria to apply and
    their limiting coefficients. It ultimately sets runTime.deltaT().

    If 'adjustTimeStep' is false, this class performs no action.

SourceFiles
    plasmaTimeControl.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaTimeControl_H
#define plasmaTimeControl_H

#include "fvCFD.H"
#include "IOdictionary.H"

#include "SoPLASMAConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class plasmaTransport;    
/*---------------------------------------------------------------------------*\
                      Class plasmaTimeControl Declaration
\*---------------------------------------------------------------------------*/

class plasmaTimeControl
{
private:

    // Private data

        //- Reference to the Time object
        Time& runTime_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary to hold control parameters
        dictionary dict_;

        // Control switches and limits

            //- Switch to enable/disable time step adjustment
            bool adjustTimeStep_;

            //- Maximum physical time step allowed
            scalar maxDeltaT_;

        // Plasma specific criteria

            //- Switch to limit Dielectric Relaxation Ratio
            Switch limitDielectricRelaxationRatio_;

            //- Switch to print Dielectric Relaxation Ratio
            Switch printDielectricRelaxationRatio_;

            //- Coefficient for Dielectric Relaxation ratio 
            scalar maxDielectricRelaxationRatio_;

            //- Switch to enable/disable Species Courant limiting
            Switch limitSpeciesCo_;

            //- Switch to print Species Courant limiting
            Switch printSpeciesCo_;

            //- Maximum allowed convective Courant for the tracked species
            scalar maxSpeciesConvectiveCo_;

            //- Maximum allowed diffusive Courant for the tracked species
            scalar maxSpeciesDiffusiveCo_;

            //- Name of the species to track for Courant limiting (default: "e")
            word courantSpeciesName_;

            //- Switch to limit by Chemistry (Source term)
            Switch limitChemistryCo_;

            //- Switch to print Chemistry source restriction
            Switch printChemistryCo_;

            //- Maximum allowed relative change due to chemistry
            scalar maxChemistryCo_;

            //- Switch to limit voltage rise rate
            Switch limitVoltageRiseRate_;

            //- Switch to print voltage rise rate info
            Switch printVoltageRiseRate_;

            //- Maximum allowed voltage rise rate [V/s]
            scalar maxVoltageRiseRate_;

            //- Name of the patch where voltage is applied
            word voltagePatchName_;

            //- Previous timestep voltage on the applied patch [V]
            scalar prevPatchVoltage_;

        // Private Member Functions

            //- Read controls from dictionary
            void read();

            //- Compute average voltage on the applied patch
            scalar patchVoltageAvg(const plasmaTransport& transport) const;

public:

    // Constructors

        //- Construct from Time and Mesh
        plasmaTimeControl(Time& runTime, const fvMesh& mesh);

    //- Destructor
    virtual ~plasmaTimeControl() {}

    // Member Functions

        //- Adjust deltaT using the specific plasmaTransport object
        //  Note: We pass it by const reference
        void adjustDeltaT(const plasmaTransport& transport);

        //- Set the first time step of the simulation
        void setInitialDeltaT(const plasmaTransport& transport);

    // Accessors (Getters)

        //- Is adjustable time stepping enabled?
        inline bool adjustTimeStep() const
        {
            return adjustTimeStep_;
        }

        //- Is dielectric relaxation limiting enabled?
        inline bool limitDielectricRelaxationRatio() const
        {
            return limitDielectricRelaxationRatio_;
        }

        //- Get the max allowed dielectric relaxation ratio
        inline scalar maxDielectricRelaxationRatio() const
        {
            return maxDielectricRelaxationRatio_;
        }

        //- Is species Courant number limiting enabled?
        inline bool limitSpeciesCo() const
        {
            return limitSpeciesCo_;
        }

        //- Get the max allowed convective species Courant number
        inline scalar maxSpeciesConvectiveCo() const
        {
            return maxSpeciesConvectiveCo_;
        }

        //- Get the max allowed diffusive species Courant number
        inline scalar maxSpeciesDiffusiveCo() const
        {
            return maxSpeciesDiffusiveCo_;
        }

        //- Get the species name used for Courant number limiting
        inline const word& courantSpeciesName() const
        {
            return courantSpeciesName_;
        }

        //- Is chemistry limiting enabled?
        inline bool limitChemistryCo() const
        {
            return limitChemistryCo_;
        }

        //- Get the max allowed chemistry limit
        inline scalar maxChemistryCo() const
        {
            return maxChemistryCo_;
        }

        //- Is voltage rise rate limiting enabled?
        inline bool limitVoltageRiseRate() const
        {
            return limitVoltageRiseRate_;
        }

        //- Get the maximum allowed voltage rise rate [V/dt]
        inline scalar maxVoltageRiseRate() const
        {
            return maxVoltageRiseRate_;
        }

        //- Get the name of the patch where voltage is applied
        inline const word& voltagePatchName() const
        {
            return voltagePatchName_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
