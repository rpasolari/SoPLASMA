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
    singleRegionElectrostaticFoam

Description
    Single-region electrostatic potential solver.

Usage
    \b singleRegionElectrostaticFoam [OPTIONS]

    Example:
        singleRegionElectrostaticFoam -case [CASE]

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "electromagneticsModel.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Single-region electrostatic potential solver."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    // No dielectric meshes for single-region
    UPtrList<fvMesh> dielectricMeshes;

    autoPtr<electromagneticsModel> em =
        electromagneticsModel::New(mesh, dielectricMeshes);

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
