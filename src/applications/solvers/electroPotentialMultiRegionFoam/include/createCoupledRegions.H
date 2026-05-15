/*---------------------------------------------------------------------------*\
  File: createCoupledRegions.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createCoupledRegions.H in OpenFOAM v2412.
      Detects and initializes coupled plasma regions based on user controls
      and boundary conditions. Creates an fvMatrix assembly for implicit
      coupling between different domains when required.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

bool coupled = false;

autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr;

forAll(fluidRegions, i)
{
    const volScalarField& ePotential = ePotentialFluid[i];
    forAll(ePotential.boundaryField(), patchI)
    {
        if (ePotential.boundaryField()[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(dielectricRegions, i)
{
    const volScalarField& ePotential = ePotentialDielectric[i];
    forAll(ePotential.boundaryField(), patchI)
    {
        if (ePotential.boundaryField()[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(fluidRegions, i)
{
    if(coupled)
    {
        Info<< "Regions are COUPLED; "
            << "ePotential will be solved MONOLITHICALLY." << nl
            << "Assembling ePotential fvMatrixAssembly from region "
            << fluidRegions[i].name() << nl << endl;

        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                ePotentialFluid[i],
                dimensionSet(0,0,1,0,0,1,0)
            )
        );
     break;
    }
}

if (coupled && !fvMatrixAssemblyPtr)
{
    forAll(dielectricRegions, i)
    {
        Info<< "Regions are COUPLED; "
            << "ePotential will be solved MONOLITHICALLY." << nl
            << "Assembling ePotential fvMatrixAssembly from region "
            << dielectricRegions[i].name() << nl << endl;

        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                ePotentialDielectric[i],
                dimensionSet(0,0,1,0,0,1,0)
            )
        );
        break;
    }
}

if (!coupled)
{
    Info<< "Regions are NOT coupled; ePotential will be solved in "
        << "SEGREGATE way" << nl << endl;
}
