/*---------------------------------------------------------------------------*\
  File: createMeshes.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createMeshes.H in OpenFOAM v2412.
      Creates and initializes the region meshes.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

regionProperties rp(runTime);

#include "createGasMesh.H"
#include "createDielectricMeshes.H"

const label nTotalRegions = gasNames.size() + dielectricNames.size();
