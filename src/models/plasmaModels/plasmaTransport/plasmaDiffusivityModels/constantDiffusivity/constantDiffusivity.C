/*---------------------------------------------------------------------------*\
  File: constantDiffusivity.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::constantDiffusivity.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "constantDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(constantDiffusivity, 0);
addToRunTimeSelectionTable
(
    plasmaDiffusivityModel,
    constantDiffusivity,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantDiffusivity::constantDiffusivity
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaDiffusivityModel(modelName, dict, mesh, species, specieIndex),
    D0_(0.0)
{
    if (!dict_.found("D"))
    {
        FatalIOErrorInFunction(dict_)
            << "Diffusivity model '" << modelName
            << "' requires entry 'D' in diffusivityCoeffs dictionary."
            << exit(FatalIOError);
    }

    dict_.lookup("D") >> D0_;

    if (D0_ < 0)
    {
        FatalIOErrorInFunction(dict)
        << "Diffusivity 'D' must be non-negative" 
        << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void constantDiffusivity::correct(volScalarField& D) const
{

    D = dimensionedScalar
    (
        "D",
        D.dimensions(),
        D0_
    );

    D.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
