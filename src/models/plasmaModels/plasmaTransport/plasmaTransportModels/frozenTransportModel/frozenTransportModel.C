/*---------------------------------------------------------------------------*\
  File: frozenTransportModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::frozenTransportModel.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "frozenTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(frozenTransportModel, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    frozenTransportModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frozenTransportModel::frozenTransportModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    plasmaTransportModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex, 
        E
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void frozenTransportModel::correct(const surfaceScalarField& phiE)
{
    // Do nothing here    
}

tmp<fvScalarMatrix> frozenTransportModel::nEqn() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    return tmp<fvScalarMatrix>(fvm::ddt(n));
}

void frozenTransportModel::updateParticleFlux(surfaceScalarField& flux) const
{
    // flux *= 0.0;
}

tmp<surfaceScalarField> frozenTransportModel::phi() const
{
    return *zeroSurfaceFieldPtr_;
}

const volVectorField& frozenTransportModel::driftVelocity() const
{
    return *zeroVectorFieldPtr_;
}

tmp<volScalarField> frozenTransportModel::electricalConductivity() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "sigma0", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    );
}

tmp<volScalarField> frozenTransportModel::diffusiveChargeSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "diffSrc0", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0.0)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
