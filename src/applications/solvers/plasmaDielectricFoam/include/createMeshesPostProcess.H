/*---------------------------------------------------------------------------*\
  File: createMeshesPostProcess.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createMeshesPostProcess.H in OpenFOAM v2412.
      Includes mesh creation and selects a representative region mesh for
      post-processing operations. Ensures that at least one mesh exists
      before proceeding.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "createMeshes.H"

if (!dielectricRegions.size())
{
    FatalErrorIn(args.executable())
        << "No region meshes present" << exit(FatalError);
}

fvMesh& mesh = gasMesh().nCells() ? gasMesh() : dielectricRegions[0];
