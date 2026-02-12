/*---------------------------------------------------------------------------*\
  File: constantMobility.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::constantMobility.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "constantMobility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(constantMobility, 0);
addToRunTimeSelectionTable(plasmaMobilityModel, constantMobility, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantMobility::constantMobility
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaMobilityModel(modelName, dict, mesh, species, specieIndex),
    mu0_(0.0)
{
    if (!dict_.found("mu"))
    {
        FatalIOErrorInFunction(dict_)
            << "Mobility model '" << modelName
            << "' requires entry 'mu' in mobilityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("mu") >> mu0_;

    if (mu0_ < 0)
    {
        FatalIOErrorInFunction(dict)
        << "Mobility 'mu' must be non-negative" 
        << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void constantMobility::correct(volScalarField& mu) const
{

    mu = dimensionedScalar
    (
        "mu",
        mu.dimensions(),
        mu0_
    );

    mu.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
