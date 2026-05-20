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

    if (plasmaDict.found("plasmaSimulationDiagnostics"))
    {
        dict_ = plasmaDict.subDict("plasmaSimulationDiagnostics");

        printSpecies_ =
            dict_.lookupOrDefault<Switch>("printSpecies", true);

        printElectromagnetics_ =
            dict_.lookupOrDefault<Switch>("printElectromagnetics", true);
    }
}

void plasmaSimulationDiagnostics::report()
{
    read();

    if (!printSpecies_ && !printElectromagnetics_)
        return;

    auto fmtRow = [](
        const std::string& label,
        const std::string& unit,
        scalar             minVal,
        scalar             maxVal,
        int                lw = 24) -> std::string
    {
        std::string line = "  " + label + " [" + unit + "]:";
        while (static_cast<int>(line.size()) < lw) line += ' ';
        line += "min = " + std::string(Foam::name(minVal).c_str());
        while (static_cast<int>(line.size()) < lw + 24) line += ' ';
        line += "max = " + std::string(Foam::name(maxVal).c_str());
        return line;
    };

    Info<< nl << "  Plasma Diagnostics" << nl
        << "  " << std::string(52, '-').c_str() << nl;

    if (printSpecies_)
    {
        Info<< "  Species number densities" << nl;
        for (label i = 0; i < transport_.species().nSpecies(); ++i)
        {
            const volScalarField& n = transport_.species().numberDensity(i);
            const word& name = transport_.species().speciesName(i);

            Info<< "  " << fmtRow
                (
                    name.c_str(), "m^-3",
                    gMin(n), gMax(n)
                ).c_str() << nl;
        }
    }

    if (printElectromagnetics_)
    {
        
        const volScalarField& Emag = transport_.species().em().Emag();
        const volScalarField& surfCharge = 
                                         transport_.species().em().surfCharge();
        
        scalar sigmaMin =  GREAT;
        scalar sigmaMax = -GREAT;
        forAll(surfCharge.boundaryField(), patchi)
        {
            const scalarField& sp = surfCharge.boundaryField()[patchi];
            if (sp.size())
            {
                sigmaMin = min(sigmaMin, min(sp));
                sigmaMax = max(sigmaMax, max(sp));
            }
        }
        reduce(sigmaMin, minOp<scalar>());
        reduce(sigmaMax, maxOp<scalar>());
        if (sigmaMin > sigmaMax) { sigmaMin = sigmaMax = 0; }

        Info<< "  Electromagnetics" << nl;

        Info<< "  " << fmtRow
            (
                "Emag",       "V/m",
                gMin(Emag),   gMax(Emag)
            ).c_str() << nl;

        Info<< "  " << fmtRow
            (
                "surfCharge",      "C/m^2",
                sigmaMin,  sigmaMax
            ).c_str() << nl;
    }

    Info<< "  " << std::string(52, '-').c_str() << nl << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
