/*---------------------------------------------------------------------------*\
  File: plasmaMobilityModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaMobilityModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaMobilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaMobilityModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaMobilityModel::plasmaMobilityModel
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    modelType_(modelType),
    mesh_(mesh),
    species_(species),
    specieIndex_(specieIndex),
    dict_(dict),
    isUniform_(false),
    isConstant_(false),
    evaluator_
    (
        plasmaPropertyEvaluator::New
        (
            modelType_,
            dict_,
            mesh_,
            "mu_" + species.speciesName(specieIndex),
            dimensionSet(-1, 0, 2, 0, 0, 1, 0)
        )
    ),
    mu_
    (
        IOobject
        (
            "mu_" + species.speciesName(specieIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(-1, 0, 2, 0, 0, 1, 0),
            0.0
        )
    )
{
    isUniform_ = evaluator_->isUniform();
    isConstant_ = evaluator_->isConstant();

    evaluator_->correct(mu_);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaMobilityModel> plasmaMobilityModel::New
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
{
    return autoPtr<plasmaMobilityModel>
    (
        new plasmaMobilityModel(modelType, dict, mesh, species, specieIndex)
    );
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaMobilityModel::correct()
{
    if (isConstant_)
    {
        return;
    }

    evaluator_->correct(mu_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
