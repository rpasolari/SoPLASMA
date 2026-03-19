/*---------------------------------------------------------------------------*\
  File: plasmaTransportModelI.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaTransportModel class.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransportModelI_H
#define Foam_plasmaTransportModelI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return model name
inline const Foam::word& 
Foam::plasmaTransportModel::modelName() const
{
    return modelName_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
