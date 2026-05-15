/*---------------------------------------------------------------------------*\
  File: createGasMeshes.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createFluidMeshes.H in OpenFOAM v2412.
      Creates and initializes the gas-region mesh.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

autoPtr<dynamicFvMesh> gasMesh;

Info<< "Create gas mesh for region " << gasRegionName
    << " for time = " << runTime.timeName() << nl << endl;

gasMesh.reset
(
    dynamicFvMesh::New
    (
        IOobject
        (
            gasRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ).ptr()
);

// Handle Dry-Run flags separately
if (args.dryRun() || args.found("dry-run-write"))
{
    Info<< "Operating in 'dry-run' mode" << endl;
    FieldBase::allowConstructFromLargerSize = true;

    if (args.found("dry-run-write"))
    {
        runTime.stopAt(Time::saWriteNow);
        gasMesh().setInstance(runTime.timeName());
    }
    else
    {
        runTime.stopAt(Time::saNoWriteNow);
    }

    functionObject::outputPrefix = "postProcessing-dry-run";
}
