/*---------------------------------------------------------------------------*\
  File: plasmaEnergyModelI.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaEnergyModel class.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaEnergyModelI_H
#define Foam_plasmaEnergyModelI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return model name
inline const Foam::word& 
Foam::plasmaEnergyModel::modelName() const
{
    return modelName_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
