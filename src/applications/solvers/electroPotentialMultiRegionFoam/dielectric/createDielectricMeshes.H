/*---------------------------------------------------------------------------*\
  File: createDielectricMeshes.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createSolidMeshes.H in OpenFOAM v2412.
      Creates and initializes the dielectric-region meshes.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

const wordList dielectricNames(rp["dielectric"]);

PtrList<fvMesh> dielectricRegions(dielectricNames.size());

forAll(dielectricNames, i)
{
    if (args.dryRun() || args.found("dry-run-write"))
    {
        Info
            << "Operating in 'dry-run' mode:"
            << " case will run for 1 time step.  "
            << "All checks assumed OK on a clean exit" << endl;

        FieldBase::allowConstructFromLargerSize = true;

        // Create a simplified 1D mesh and attempt to re-create boundary
        // conditions
        dielectricRegions.set
        (
            i,
            new simplifiedMeshes::columnFvMesh(runTime, dielectricNames[i])
        );

        // Stop after 1 iteration of the simplified mesh
        if (args.found("dry-run-write"))
        {
            // Using saWriteNow triggers function objects execute(), write()
            runTime.stopAt(Time::saWriteNow);

            // Make sure mesh gets output to the current time (since
            // instance no longer constant)
            dielectricRegions[i].setInstance(runTime.timeName());
        }
        else
        {
            // Using saNoWriteNow triggers function objects execute(),
            // but not write()
            runTime.stopAt(Time::saNoWriteNow);
        }

        functionObject::outputPrefix = "postProcessing-dry-run";
    }    
    else
    {
        Info << "Create dielectric mesh for region " << dielectricNames[i]
             << " for time = " << runTime.timeName() << nl << endl;

        dielectricRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    dielectricNames[i],
                    runTime.timeName(),
                    runTime,
                    IOdictionary::MUST_READ
                )
            )
        );
    }
}
