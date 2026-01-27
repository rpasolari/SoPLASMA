/*---------------------------------------------------------------------------*\
  File: solveElectricPotential.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Selects and executes the electrostatic potential solver (coupled or 
      non-coupled) for the current iteration.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if(coupled)
{
    #include "solveElectricPotentialCoupled.H"
}
else
{
    #include "solveElectricPotentialNonCoupled.H"
}





