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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    electroPotentialMultiRegionFoam

Description
    Steady-state multi-region electrostatic potential solver for coupled
    fluid and dielectric domains.

Usage
    \b electroPotentialMultiRegionFoam [OPTIONS]

    Example:
        electroPotentialMultiRegionFoam -case [CASE]

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "SoPLASMAConstants.H"
#include "fvCFD.H"
#include "regionProperties.H"
#include "coordinateSystem.H"
#include "loopControl.H"

#include "fvSolution.H"
#include "solutionControl.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state multi-region electrostatic potential solver for coupled"
        " fluid and dielectric domains."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    #include "createCoupledRegions.H"

    #include "readElectricPotentialControls.H"

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;
        
        // Solve the Poisson/Laplace Equation (electric potential)
        if(coupled)
        {
            #include "solveElectricPotentialCoupled.H"
        }
        else
        {
            #include "solveElectricPotentialNonCoupled.H"
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
