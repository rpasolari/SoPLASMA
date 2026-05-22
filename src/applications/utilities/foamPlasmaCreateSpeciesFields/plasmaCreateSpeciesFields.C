/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    plasmaCreateSpeciesFields

Description
    Utility to generate template number-density fields (n_speciesName) in the
    0/ directory for all species listed in the plasmaSpecies dictionary.

    The tool:
      - detects and uses the gas region (if present in a multi-region case),
      - reads species from the plasmaSpecies dictionary,
      - loads the mesh of the region,
      - creates volScalarField files with default values and zeroGradient
        boundary conditions for all patches.

    These generated fields are **only templates**.
    The user must review and modify the internalField values, boundary
    conditions, and any other field settings according to the needs of the
    specific plasma simulation.

Usage
    \b plasmaDielectricFoam [OPTIONS]

    Example:
        plasmaCreateSpeciesFields -case plasmaCase

        or 

        mpirun -np 4 plasmaCreateFields -parallel

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Create species fields in 0 folder"
    );

    argList::noBanner();
    argList::noJobInfo();
    argList::noFunctionObjects();

    argList::addOption
    (
        "dict",
        "file",
        "Alternative plasmaSpecies dictionary"
    );

    argList args(argc, argv);

    // Create minimal Time
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    word gasRegionName;

    fileName globalConstantDir = runTime.globalPath()/"constant";
    fileName regionPropsFile = runTime.constant()/"regionProperties";

    // Find the gas region, or default region in single region cases
    if (isFile(regionPropsFile))
    {
        if (Pstream::master())
        {
            Info << "Found regionProperties: " << regionPropsFile << nl;
        }

        regionProperties rp(runTime);

        if (rp.found("gas"))
        {
            const wordList& names = rp["gas"];

            if (names.size() == 0)
            {
                FatalErrorInFunction
                    << "A 'gas' region type is defined in regionProperties, "
                    << "but no region names are listed under it." << nl
                    << exit(FatalError);
            }

            if (names.size() > 1)
            {
                FatalErrorInFunction
                    << "Multiple gas regions detected in "
                    << "constant/regionProperties:" << nl
                    << "  gas: " << names << nl
                    << "This utility supports ONLY ONE gas region." << nl
                    << exit(FatalError);
            }

            // Exactly one gas region
            gasRegionName = names[0];

            if (Pstream::master())
            {
                Info << "Detected gas region: " << gasRegionName << nl;
            }
        }
    }

    // Report region used
    if (Pstream::master())
    {
        if (gasRegionName.empty())
            Info << "No gas region found → using constant/ as path." << nl;
        else
            Info << "Using region: " << gasRegionName << nl;
    }

    fileName plasmaSpeciesPath;

    // Find the plasmaSpecies file in the case
    if (!gasRegionName.empty())
    {
        // Multi-region path
        fileName regionDir = globalConstantDir/gasRegionName;

        if (isFile(regionDir/"plasmaSpeciesProperties"))
        {
            plasmaSpeciesPath = regionDir/"plasmaSpeciesProperties";
        }
        else
        {
            FatalErrorInFunction
                << "Region '" << gasRegionName << "' exists but has no "
                << "plasmaSpeciesProperties file:" << nl
                << "  " << regionDir/"plasmaSpeciesProperties" << nl
                << exit(FatalError);
        }
    }
    else
    {
        // Single-region path
        if (isFile(runTime.constant()/"plasmaSpeciesProperties"))
        {
            plasmaSpeciesPath = runTime.constant()/"plasmaSpeciesProperties";
        }
        else
        {
            FatalErrorInFunction
                << "No 'constant/plasmaSpeciesProperties' file found." << nl 
                << exit(FatalError);
        }
    }

    // Read plasmaSpecies dict
    if (Pstream::master())
    {
        Info << "Reading plasmaSpecies from: " << plasmaSpeciesPath << nl;
    }

    IOdictionary speciesDict
    (
        IOobject
        (
            "plasmaSpeciesProperties",
            plasmaSpeciesPath.path(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!speciesDict.found("activeSpecies"))
    {
        FatalErrorInFunction
            << "The plasmaSpecies dictionary has no 'activeSpecies' entry!"
            << exit(FatalError);
    }

    wordList species(speciesDict.lookup("activeSpecies"));

    if (Pstream::master())
    {
        Info << "Found " << species.size() << " active species." << nl;
    }

    for (const word& s : species)
    {
        Info << "  - " << s << nl;
    }

    // Determine region path for mesh creation
    word regionToLoad = gasRegionName;

    // Create the mesh
    fvMesh mesh
    (
        IOobject
        (
            regionToLoad,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Field Creation Loop
    dimensionSet dimDensity(0, -3, 0, 0, 0, 0, 0);

    for (const word& s : species)
    {
        const word fieldName = "n_" + s;

        if (Pstream::master())
        {
            Info << "Creating field: " << fieldName << endl;
        }

        // Create the field using the specific constructor that sets:
        // 1. Internal Field Value (0)
        // 2. Physical Boundary Type (zeroGradient)
        // 3. Processor Boundary Type (calculated automatically)
        volScalarField nSpecies
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity, 0.0),
            word("zeroGradient")
        );

        // Write the field
        // In serial: writes to case/0/
        // In parallel: writes to case/processorN/0/
        nSpecies.write();
    }

    if (Pstream::master())
    {
        Info << nl << "Fields created successfully." << nl;
    }

    return 0;
}
