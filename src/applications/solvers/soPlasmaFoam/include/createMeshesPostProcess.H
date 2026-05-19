/*---------------------------------------------------------------------------*\
  File: createMeshesPostProcess.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createMeshesPostProcess.H in OpenFOAM v2412.
      Includes mesh creation and selects a representative region mesh for
      post-processing operations.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "createMeshes.H"

if (!gasMesh().nCells() && !dielectricRegions.size())
{
    FatalErrorIn(args.executable())
        << "This processor has no mesh cells in any region"
        << exit(FatalError);
}

fvMesh& mesh = gasMesh().nCells() ? gasMesh() : dielectricRegions[0];
