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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    soPlasmaFoam

Description
    Transient solver for coupled gas (plasma) and dielectric domains developed
    for plasma simulation purposes.

Usage
    \b soPlasmaFoam [OPTIONS]

    Example:
        soPlasmaFoam -case testCase

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include <iomanip>

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "loopControl.H"
#include "fvSolution.H"
#include "solutionControl.H"
#include "mappedPatchBase.H"
#include "pimpleControl.H"

#include "plasmaConstants.H"
#include "plasmaTimeControl.H"
#include "electromagneticsModel.H"
#include "multiRegionPoissonModel.H"
#include "dielectricRegion.H"
#include "plasmaSpecies.H"
#include "plasmaTransport.H"
#include "driftDiffusion.H"
#include "plasmaSimulationDiagnostics.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for coupled gas (plasma) and dielectric domains"
        " developed for plasma simulation purposes."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    //- Create the electromagnetics model
    autoPtr<electromagneticsModel> em =
        electromagneticsModel::New(gasMesh(), dielectricFvMeshes);

    //- Create the plasmaSpecies model
    plasmaSpecies species(gasMesh(), em());

    //- Create the plasmaTransport model
    plasmaTransport transport(gasMesh(), species);

    //- Create the PIMPLE loop control
    pimpleControl pimple(gasMesh());

    //- Create the timeControl manager and set initial time-step
    plasmaTimeControl timeControl(runTime, gasMesh());
    timeControl.setInitialDeltaT(transport);

    //- Create the plasmaSimulationDiagnostics manager
    plasmaSimulationDiagnostics diagnostics(runTime, transport);

    #include "reportSimulationSummary.H"

    runTime.writeNow();

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        timeControl.adjustDeltaT(transport);

        // Semi-implicit Poisson branch
        if (em->PoissonScheme() == "semiImplicit")
        {
            while (pimple.loop())
            {
                // Solve electromagnetics
                em->solve
                (
                    transport.electricalConductivity(), 
                    transport.diffusiveChargeSource()
                );

                // Solve transport equations
                transport.solve();

                if (pimple.finalIter())
                {
                    // Update charge density
                    species.updateChargeDensity();

                    // Update surface charge
                    transport.updateSurfaceCharge();
                }
            }
        }
        else //Explicit Poisson branch
        {
            while (pimple.loop())
            {
                // Solve transport equations
                transport.solve();

                // Update charge density
                species.updateChargeDensity();

                // Update surface charge
                transport.updateSurfaceCharge();

                em->solve();
            }
        }
        
        diagnostics.report();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
