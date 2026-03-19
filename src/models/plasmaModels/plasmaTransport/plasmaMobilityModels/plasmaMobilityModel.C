/*---------------------------------------------------------------------------*\
  File: plasmaMobilityModel.C
  Part of: foamPlasmaToolkit
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

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void plasmaMobilityModel::update() const
{
    // Skip recalculation if the model is constant in time and already set
    if (isConstant_ && upToDate_) return;

    if (isUniform_)
    {
        // Update the single scalar value
        muValue_ = evaluator_->evaluateScalar();
    }
    else
    {
        // Create or update the volScalarField
        if (!muFieldPtr_.valid())
        {
            muFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject("mu", mesh_.time().timeName(), mesh_),
                    mesh_, muValue_.dimensions()
                )
            );
        }
        evaluator_->correct(muFieldPtr_.ref());
    }

    upToDate_ = true;
}

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
    upToDate_(false),
    muValue_("mu", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
{

    evaluator_ = plasmaPropertyEvaluator::New
    (
        modelType_,
        dict_,
        mesh_,
        species_,
        specieIndex_,
        "Mobility",
        dimensionSet(-1, 0, 2, 0, 0, 1, 0)
    );

    isUniform_ = evaluator_->isUniform();
    isConstant_ = evaluator_->isConstant();
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

tmp<volScalarField> plasmaMobilityModel::mu() const
{
    update();

    if (isUniform_)
    {
        // Return a virtual field
        return tmp<volScalarField>::New
        (
            IOobject("mu_tmp", mesh_.time().timeName(), mesh_),
            mesh_, muValue_
        );
    }

    // Return reference to the cached field
    return *muFieldPtr_;
}

void plasmaMobilityModel::correct(volScalarField& mu) const
{
    update();

    if (isUniform_) mu = muValue_;
    else mu = *muFieldPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
