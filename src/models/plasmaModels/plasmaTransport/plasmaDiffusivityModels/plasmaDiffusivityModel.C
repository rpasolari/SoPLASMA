/*---------------------------------------------------------------------------*\
  File: plasmaDiffusivityModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaDiffusivityModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaDiffusivityModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaDiffusivityModel::plasmaDiffusivityModel
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
            "D_" + species.speciesName(specieIndex),
            dimensionSet(0, 2, -1, 0, 0, 0, 0)
        )
    ),
    D_
    (
        IOobject
        (
            "D_" + species.speciesName(specieIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            0.0
        )
    )
{
    isUniform_ = evaluator_->isUniform();
    isConstant_ = evaluator_->isConstant();

    evaluator_->correct(D_);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaDiffusivityModel> plasmaDiffusivityModel::New
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
{
    return autoPtr<plasmaDiffusivityModel>
    (
        new plasmaDiffusivityModel(modelType, dict, mesh, species, specieIndex)
    );
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaDiffusivityModel::correct()
{
    if (isConstant_)
    {
        return;
    }

    evaluator_->correct(D_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
