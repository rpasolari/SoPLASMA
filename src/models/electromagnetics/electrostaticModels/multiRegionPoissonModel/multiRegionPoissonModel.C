/*---------------------------------------------------------------------------*\
  File: singleRegionPoissonModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::singleRegionPoissonModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "singleRegionPoissonModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(singleRegionPoissonModel, 0);
addToRunTimeSelectionTable
(
    electromagneticsModel,
    singleRegionPoissonModel,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// CONSTRUCTOR - No "..." here!
multiRegionPoissonModel::multiRegionPoissonModel
(
    const dynamicFvMesh& mesh,
    const UPtrList<dynamicFvMesh>& dielectricMeshes,
    const plasmaSpecies& species,
    const plasmaTransport& transport,
    const dictionary& dict
)
:
    electromagneticsModel(mesh, dielectricMeshes, species, transport, dict),
    
    ePotential_
    (
        IOobject("ePotential", mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    ),
    E_
    (
        IOobject("E", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)
    ),
    ePotentialDielectricList_(dielectricMeshes.size()),
    EDielectricList_(dielectricMeshes.size()),
    epsilonDielectricList_(dielectricMeshes.size()),
    coupled_(false),
    fvMatrixAssemblyPtr_(nullptr)
{
    // 1. Initialize Epsilon for Gas
    scalar epsilonR = dict.getOrDefault<scalar>("dielectricConstant", 1.0);
    epsilon_ = dimensionedScalar("epsilon", epsilonR * constant::plasma::epsilon0);

    // 2. Initialize Dielectric Fields and Epsilons
    forAll(dielectricMeshes_, i)
    {
        const dynamicFvMesh& dMesh = dielectricMeshes_[i];
        
        ePotentialDielectricList_.set(i, new volScalarField(
            IOobject("ePotential", dMesh.time().timeName(), dMesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), dMesh
        ));

        EDielectricList_.set(i, new volVectorField(
            IOobject("E", dMesh.time().timeName(), dMesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            dMesh, dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)
        ));

        // Assuming dielectrics have their own epsilon in the dictionary
        epsilonDielectricList_.set(i, new dimensionedScalar("epsilon", 
            constant::plasma::epsilon0)); // Update this logic to read from dict if needed
    }

    // 3. INTEGRATE YOUR COUPLING LOGIC
    // Check Gas
    forAll(ePotential_.boundaryField(), patchI)
    {
        if (ePotential_.boundaryField()[patchI].useImplicit())
        {
            coupled_ = true; break;
        }
    }

    // Check Dielectrics
    if (!coupled_)
    {
        forAll(ePotentialDielectricList_, i)
        {
            forAll(ePotentialDielectricList_[i].boundaryField(), patchI)
            {
                if (ePotentialDielectricList_[i].boundaryField()[patchI].useImplicit())
                {
                    coupled_ = true; break;
                }
            }
            if (coupled_) break;
        }
    }

    // 4. PERSISTENT MATRIX ALLOCATION (Your Logic)
    if (coupled_)
    {
        Info << "Regions are COUPLED; Initializing persistent Matrix Assembly." << endl;
        fvMatrixAssemblyPtr_.reset
        (
            new fvMatrix<scalar>
            (
                ePotential_,
                dimensionSet(0,0,1,0,0,1,0)
            )
        );
    }
}

void multiRegionPoissonModel::solve()
{
    const volScalarField& rhoq = mesh_.lookupObject<volScalarField>("chargeDensity");

    if (coupled_ && fvMatrixAssemblyPtr_)
    {
        // Re-use the persistent matrix (overwrite coefficients)
        fvMatrixAssemblyPtr_() = (fvm::laplacian(epsilon_, ePotential_) == -rhoq);

        forAll(dielectricMeshes_, i)
        {
            fvMatrixAssemblyPtr_() += 
                fvm::laplacian(epsilonDielectricList_[i], ePotentialDielectricList_[i]);
        }

        fvMatrixAssemblyPtr_->solve();
    }
    else
    {
        // Segregated solve
        fvScalarMatrix gasEqn(fvm::laplacian(epsilon_, ePotential_) == -rhoq);
        gasEqn.solve();

        forAll(dielectricMeshes_, i)
        {
            fvScalarMatrix dielEqn(fvm::laplacian(epsilonDielectricList_[i], ePotentialDielectricList_[i]));
            dielEqn.solve();
        }
    }

    // Final Field Updates
    E_ = -fvc::grad(ePotential_);
    E_.correctBoundaryConditions();
    forAll(EDielectricList_, i)
    {
        EDielectricList_[i] = -fvc::grad(ePotentialDielectricList_[i]);
        EDielectricList_[i].correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
