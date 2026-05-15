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

PtrList<dynamicFvMesh> dielectricRegions(dielectricNames.size());

forAll(dielectricNames, i)
{
    Info << "Create dielectric mesh for region " << dielectricNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

    dielectricRegions.set
    (
        i,
        dynamicFvMesh::New
        (
            IOobject
            (
                dielectricNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).ptr()
    );

    // Handle Dry-Run flags matching the gasMesh implementation
    if (args.dryRun() || args.found("dry-run-write"))
    {
        Info<< "Operating in 'dry-run' mode for region " 
            << dielectricNames[i] << endl;
        
        FieldBase::allowConstructFromLargerSize = true;

        if (args.found("dry-run-write"))
        {
            // Using saWriteNow triggers function objects execute(), write()
            runTime.stopAt(Time::saWriteNow);

            // Ensure mesh gets output to the current time step
            dielectricRegions[i].setInstance(runTime.timeName());
        }
        else
        {
            // Using saNoWriteNow triggers function objects execute(), but not write()
            runTime.stopAt(Time::saNoWriteNow);
        }

        functionObject::outputPrefix = "postProcessing-dry-run";
    }
}
