/*---------------------------------------------------------------------------*\
  File: plasmaSimulationDiagnostics.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSimulationDiagnostics.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaSimulationDiagnostics.H"
#include "plasmaTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaSimulationDiagnostics::plasmaSimulationDiagnostics
(
    const Time& runTime,
    const plasmaTransport& transport
)
:
    runTime_(runTime),
    transport_(transport),
    dict_(dictionary::null),
    printSpecies_(true),
    printElectromagnetics_(true)
{
    read();

    Info << "Plasma simulation diagnostics succesfully constructed." << nl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasmaSimulationDiagnostics::read()
{
    IOdictionary plasmaDict
        (
            IOobject
            (
                "plasmaSimulationControls",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    if (plasmaDict.found("plasmaDiagnostics"))
    {
        dict_ = plasmaDict.subDict("plasmaDiagnostics");

        printSpecies_ =
            dict_.lookupOrDefault<Switch>("printSpecies", true);

        printElectromagnetics_ =
            dict_.lookupOrDefault<Switch>("printElectromagnetics", true);
    }
}

void plasmaSimulationDiagnostics::report()
{
    // Update variables if the dictionary was modified during runtime
    read();

    // If both are false, skip printing entirely to save console space
    if (!printSpecies_ && !printElectromagnetics_)
    {
        return;
    }

    Info<< "  Plasma Diagnostics" << nl
        << "  " << std::string(52, '-').c_str() << nl;

    if (printSpecies_)
    {
        for (label i = 0; i < transport_.species().nSpecies(); ++i)
        {
            const volScalarField& n = transport_.species().numberDensity(i);
            const word& name = transport_.species().speciesName(i);
            
            Info<< "    " << name << " [m^-3] min: " << gMin(n) 
                << "  max: " << gMax(n) << nl;
        }
    }

    if (printElectromagnetics_)
    {
        const volScalarField& Emag = transport_.species().em().Emag();
        const volScalarField& surfCharge = transport_.species().em().surfCharge();

        Info<< "    Emag       [V/m] min: " << gMin(Emag) 
            << "  max: " << gMax(Emag) << nl;
            
        Info<< "    surfCharge [C/m^2] min: " << gMin(surfCharge) 
            << "  max: " << gMax(surfCharge) << nl;
    }

    Info<< "  " << std::string(52, '-').c_str() << nl << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
