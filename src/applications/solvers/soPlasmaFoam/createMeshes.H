/*---------------------------------------------------------------------------*\
  File: createMeshes.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createMeshes.H in OpenFOAM v2412.
      Creates and initializes the region meshes.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

regionProperties rp(runTime, IOobjectOption::READ_IF_PRESENT);

wordList gasNames;
if (rp.found("gas"))
{
    gasNames = rp["gas"];
}

if (gasNames.size() > 1)
{
    FatalErrorInFunction
        << "Only one gas region is supported by this solver." << nl
        << "Found " << gasNames.size() << " regions: " << gasNames << nl
        << exit(FatalError);
}

word gasRegionName;

if (gasNames.empty())
{
    gasRegionName = polyMesh::defaultRegion;
    Info<< "No gas region specified in constant/regionProperties." << nl
        << "Falling back to default mesh: " << gasRegionName << nl << endl;
}
else
{
    gasRegionName = gasNames[0];
}

wordList dielectricNames;
if (rp.found("dielectric"))
{
    dielectricNames = rp["dielectric"];
}

// const label nTotalRegions = 
//     (gasNames.empty() ? 1 : gasNames.size()) + dielectricNames.size();

#include "createGasMesh.H"
#include "createDielectricMeshes.H"
