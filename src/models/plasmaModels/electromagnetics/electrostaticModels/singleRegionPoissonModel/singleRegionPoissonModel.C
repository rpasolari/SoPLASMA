#include "singleRegionPoissonModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// 1. Register the model
defineTypeNameAndDebug(singleRegionPoissonModel, 0);
addToRunTimeSelectionTable
(
    electromagneticsModel,
    singleRegionPoissonModel,
    dictionary
);

// 2. The Constructor - MUST have the full argument list
singleRegionPoissonModel::singleRegionPoissonModel
(
    const dynamicFvMesh& mesh,
    const UPtrList<dynamicFvMesh>& dielectricMeshes,
    const plasmaSpecies& species,
    const plasmaTransport& transport,
    const dictionary& dict
)
:
    // Pass the arguments up to the parent class
    electromagneticsModel(mesh, dielectricMeshes, species, transport, dict),
    
    // Initialize your fields
    ePotential_
    (
        IOobject
        (
            "ePotential",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)
    )
{
    // Initialize your epsilon
    scalar epsilonR = dict.getOrDefault<scalar>("dielectricConstant", 1.0);
    epsilon_ = dimensionedScalar
    (
        "epsilon",
        epsilonR * constant::plasma::epsilon0
    );

    Info<< "    Poisson Solver: singleRegionPoissonModel initialized" << endl;
}

// 3. Getter for E
const volVectorField& singleRegionPoissonModel::E() const
{
    return E_;
}

// 4. Solve logic
void singleRegionPoissonModel::solve()
{
    // Your solver logic here (fvm::laplacian, etc.)
}

} // End namespace Foam
