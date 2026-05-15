/*---------------------------------------------------------------------------*\
  File: electromagneticsModelI.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the electromagneticsModel class.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_electromagneticsModelI_H
#define Foam_electromagneticsModelI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh&
Foam::electromagneticsModel::mesh() const
{
    return mesh_;
}

inline const Foam::UPtrList<Foam::fvMesh>& 
Foam::electromagneticsModel::dielectricMeshes() const
{
    return dielectricMeshes_;
}

inline const Foam::dictionary&
Foam::electromagneticsModel::dict() const
{ 
    return *this; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
