/*---------------------------------------------------------------------------*\
  File: plasmaDiffusivityModel.C
  Part of: foamPlasmaToolkit
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

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void plasmaDiffusivityModel::update() const
{
    // Skip recalculation if the model is constant in time and already set
    if (isConstant_ && upToDate_) return;

    if (isUniform_)
    {
        // Update the single scalar value
        DValue_ = evaluator_->evaluateScalar();
    }
    else
    {
        // Create or update the volScalarField
        if (!DFieldPtr_.valid())
        {
            DFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject("D", mesh_.time().timeName(), mesh_),
                    mesh_, DValue_.dimensions()
                )
            );
        }
        evaluator_->correct(DFieldPtr_.ref());
    }

    upToDate_ = true;
}

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
    upToDate_(false),
    DValue_("D", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
{
    evaluator_ = plasmaPropertyEvaluator::New
    (
        modelType_,
        dict_,
        mesh_,
        species_,
        specieIndex_,
        "Diffusivity",
        dimensionSet(0, 2, -1, 0, 0, 0, 0)
    );

    isUniform_ = evaluator_->isUniform();
    isConstant_ = evaluator_->isConstant();
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

tmp<volScalarField> plasmaDiffusivityModel::D() const
{
    update();

    if (isUniform_)
    {
        // Return a virtual field
        return tmp<volScalarField>::New
        (
            IOobject("D_tmp", mesh_.time().timeName(), mesh_),
            mesh_, DValue_
        );
    }

    // Return reference to the cached field (Case 2 & 4)
    return *DFieldPtr_;
}

void plasmaDiffusivityModel::correct(volScalarField& D) const
{
    update();

    if (isUniform_) D = DValue_;
    else D = *DFieldPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
