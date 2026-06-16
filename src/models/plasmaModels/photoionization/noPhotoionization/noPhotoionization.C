/*---------------------------------------------------------------------------*\
  File: noPhotoionization.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::noPhotoionization.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "noPhotoionization.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(noPhotoionization, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    noPhotoionization,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPhotoionization::noPhotoionization(const fvMesh& mesh)
:
    photoionizationModel(mesh)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noPhotoionization::correct()
{
    // Do nothing
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
