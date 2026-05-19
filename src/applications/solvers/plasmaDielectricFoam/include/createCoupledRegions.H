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

if (nTotalRegions == 1)
{
    Info<< nl << "Single region simulation. "
        << "Skipping region coupling." << endl;

    coupled = false;
}
else
{
    Info<< nl << "Multiple regions detected. "
        << "Initializing coupled interfaces." << endl;

    forAll(ePotential.boundaryField(), patchI)
    {
        if (ePotential.boundaryField()[patchI].useImplicit())
        {
            coupled = true;
            break;
        }
    }

    if (!coupled)
    {
        forAll(dielectricRegions, i)
        {
            const volScalarField& ePotentialDielectric = 
                ePotentialDielectricList[i];

            forAll(ePotentialDielectric.boundaryField(), patchI)
            {
                if (ePotentialDielectric.boundaryField()[patchI].useImplicit())
                {
                    coupled = true;
                    break;
                }
            }
            if (coupled) break;
        }
    }

    if (coupled)
    {
        Info<< "Regions are COUPLED; "
            << "ePotential will be solved MONOLITHICALLY." << nl
            << "Assembling ePotential fvMatrixAssembly from region "
            << gasRegionName << nl << endl;

        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                ePotential,
                dimensionSet(0,0,1,0,0,1,0)
            )
        );
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
                    ePotentialDielectricList[i],
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
}




