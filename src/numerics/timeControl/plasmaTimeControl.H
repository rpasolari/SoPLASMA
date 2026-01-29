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
    Foam::plasmaTimeControl

Description
    Class for automatic time-step control in plasma simulations.

    This class serves as a specialized replacement for the standard OpenFOAM
    'setDeltaT.H' functionality. It calculates the maximum permissible time
    step based on plasma-specific stability criteria.

    The class reads the 'plasmaTimeControl' dictionary (typically found within
    system/controlDict) to determine which criteria to apply and their 
    limiting coefficients. It ultimately sets runTime.deltaT().

    If 'adjustTimeStep' is false, this class performs no action.

SourceFiles
    plasmaTimeControl.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaTimeControl_H
#define plasmaTimeControl_H

#include "fvCFD.H"
#include "IOdictionary.H"

#include "foamPlasmaToolkitConstants.H"

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

            //- Maximum allowed Courant number for the species
            scalar maxSpeciesCo_;

            //- Name of the species to track (default: "e")
            word speciesName_;

        // Private Member Functions

            //- Read controls from dictionary
            void read();

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

        //- Get the max allowed species Courant number
        inline scalar maxSpeciesCo() const
        {
            return maxSpeciesCo_;
        }

        //- Get the species name used for Courant number limiting
        inline const word& speciesName() const
        {
            return speciesName_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
