/*---------------------------------------------------------------------------*\
  File: plasmaTransportModelI.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaTransportModel class.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransportModelI_H
#define Foam_plasmaTransportModelI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word&
Foam::plasmaTransportModel::modelName() const
{
    return modelName_;
}

inline const Foam::fvMesh&
Foam::plasmaTransportModel::mesh() const
{
    return mesh_;
}

inline const Foam::plasmaSpecies&
Foam::plasmaTransportModel::species() const
{
    return species_;
}

inline Foam::label
Foam::plasmaTransportModel::specieIndex() const
{
    return specieIndex_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
