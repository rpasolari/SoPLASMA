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
    multiRegionElectrostaticFoam

Description
    Multi-region electrostatic potential solver.

Usage
    \b multiRegionElectrostaticFoam [OPTIONS]

    Example:
        multiRegionElectrostaticFoam -case [CASE]

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "dynamicFvMesh.H"

#include "electromagneticsModel.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Multi-region electrostatic potential solver."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"

    UPtrList<fvMesh> dielectricMeshes(dielectricRegions.size());

    forAll(dielectricRegions, i)
    {
        dielectricMeshes.set(i, &dielectricRegions[i]);
    }

    autoPtr<electromagneticsModel> em =
        electromagneticsModel::New(gasMesh(), dielectricMeshes);

    Info<< "\nCalculating electric potential distribution\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;
        
        em->solve();
    
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
