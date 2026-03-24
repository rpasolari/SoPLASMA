/*---------------------------------------------------------------------------*\
  File: multiRegionPoissonModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::multiRegionPoissonModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "mappedPatchBase.H"

#include "multiRegionPoissonModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(multiRegionPoissonModel, 0);
addToRunTimeSelectionTable
(
    electromagneticsModel,
    multiRegionPoissonModel,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
multiRegionPoissonModel::multiRegionPoissonModel
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes,
    const plasmaSpecies* species
)
:
    electromagneticsModel(mesh, dielectricMeshes, species),
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
        dimensionedVector
        (
            "zero", 
            dimensionSet(1, 1, -3, 0, 0, -1, 0), 
            vector::zero
        )
    ),
    Emag_
    (
        IOobject
        (
            "Emag", 
            mesh.time().timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mag(E_)
    ),
    phiE_
    (
        IOobject
        (
            "phiE", 
            mesh.time().timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        -fvc::snGrad(ePotential_) * mesh.magSf()
    ),
    epsilon_
    (
        "epsilon", 
        dimensionSet(-1, -3, 4, 0, 0, 2, 0), 
        0.0
    ),
    epsilonR_(1.0),
    PoissonScheme_("explicit"),
    transportPtr_(nullptr),
    dielectrics_(dielectricMeshes.size()),
    coupled_(false),
    fvMatrixAssemblyPtr_(nullptr)
{
    const word coeffsName(type() + "Coeffs");

    if (!found(coeffsName))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary '" << coeffsName << "' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& coeffs(subDict(coeffsName));
    epsilonR_ = coeffs.get<scalar>("dielectricConstant");
    epsilon_ = epsilonR_ * constant::plasma::epsilon0;

    PoissonScheme_ = coeffs.get<word>("PoissonScheme");

    const dictionary& nonCoupledResidualControl = 
        coeffs.subOrEmptyDict("nonCoupledResidualControl");

    maxNonCoupledIterations_ = 
        nonCoupledResidualControl.getOrDefault<label>("maxIter", 100);

    nonCoupledTolerance_ = 
        nonCoupledResidualControl.getOrDefault<scalar>("tolerance", 1e-6);

    if (PoissonScheme_ != "explicit" && PoissonScheme_ != "semiImplicit")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown PoissonScheme '" << PoissonScheme_ << "'." << nl
            << "Valid options are: (explicit | semiImplicit)" << nl
            << exit(FatalIOError);
    }

    forAll(dielectricMeshes_, i)
    {
        const word& regionName = dielectricMeshes_[i].name();
        if (!coeffsDict.found(regionName))
        {
            FatalIOErrorInFunction(coeffsDict)
                << "Region '" << regionName << "' not found in "
                << coeffsName << " dictionary." << nl 
                << "Every dielectric mesh must have a corresponding "
                << "sub-dictionary." << exit(FatalIOError);
        }
        const dictionary& regionDict = coeffsDict.subDict(regionName);

        dielectrics_.set
        (
            i,
            new dielectricRegion(dielectricMeshes_[i], regionDict)
        );
    }   

    bool anyImplicit = false;

    forAll(ePotential_.boundaryField(), patchI)
    {
        if (ePotential_.boundaryField()[patchI].useImplicit())
        {
            anyImplicit = true;
            break;
        }
    }

    if (!anyImplicit)
    {
        for (const auto& reg : dielectrics_)
        {
            forAll(reg.ePotential().boundaryField(), patchI)
            {
                if (reg.ePotential().boundaryField()[patchI].useImplicit())
                {
                    anyImplicit = true;
                    break;
                }
            }
            if (anyImplicit) break;
        }
    }

    if (anyImplicit)
    {
        coupled_ = true;
        Info<< "    Implicit coupling detected. "
            << "Validating interface consistency..." << endl;

        forAll(ePotential_.boundaryField(), patchI)
        {
            const fvPatch& p = ePotential_.boundaryField()[patchI].patch();

            if (isA<mappedPatchBase>(p)
             && !ePotential_.boundaryField()[patchI].useImplicit())
            {
                FatalErrorInFunction
                    << "Mixed coupling detected on patch '" << p.name() 
                    << "'." << nl << "In Monolithic mode, ALL interface "
                    << "patches must set 'useImplicit true'." << nl
                    << exit(FatalError);
            }
        }

        for (const auto& reg : dielectrics_)
        {
            forAll(reg.ePotential().boundaryField(), patchI)
            {
                const fvPatch& p = 
                    reg.ePotential().boundaryField()[patchI].patch();

                if (isA<mappedPatchBase>(p)
                 && !reg.ePotential().boundaryField()[patchI].useImplicit())
                {
                    FatalErrorInFunction
                        << "Mixed coupling detected in region '" 
                        << reg.mesh().name() << "' on patch '" << p.name() 
                        << "'." << nl << "All mapped patches must be "
                        << "consistent (all implicit or all explicit)."
                        << exit(FatalError);
                }
            }
        }

        Info<< "Regions are COUPLED; "
            << "ePotential will be solved MONOLITHICALLY." << nl
            << "Assembling ePotential fvMatrixAssembly" << endl;

        fvMatrixAssemblyPtr_.reset
        (
            new fvMatrix<scalar>(ePotential_, ePotential_.dimensions())
        );
    }
    else
    {
        coupled_ = false;
        Info<< "Regions are NOT coupled; ePotential will be solved in "
            << "SEGREGATE way" << nl << endl;
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void multiRegionPoissonModel::solve()
{
    if (PoissonScheme_ == "semiImplicit" && transportPtr_ && hasPlasma())
    {
        this->solve(*transportPtr_);
        return;
    }

    if (coupled_)
    {
        fvScalarMatrix ePotentialGasEqn
        (
            fvm::laplacian(epsilon_, ePotential_)
        );

        if (hasPlasma())
        {
            ePotentialGasEqn == -species().chargeDensity();
        }

        fvMatrixAssemblyPtr_->addFvMatrix(gasEqn);

        for (dielectricRegion& reg : dielectrics_)
        {
            fvScalarMatrix ePotentialDielectricEqn
            (
                fvm::laplacian(reg.epsilon(), reg.ePotential())
            );
            fvMatrixAssemblyPtr_->addFvMatrix(ePotentialDielectricEqn);
        }

        fvMatrixAssemblyPtr_->solve();

        ePotential_.correctBoundaryConditions();
        for (dielectricRegion& reg : dielectrics_)
        {
            reg.ePotential().correctBoundaryConditions();
        }

        fvMatrixAssemblyPtr_->clear();
    }
    else
    {
        for (int iter = 1; iter <= maxNonCoupledIterations_; ++iter)
        {
            if (maxNonCoupledIterations_ > 1)
            {
                Info<< "  -ePotential iteration (" << iter << "/" 
                    << maxNonCoupledIterations_ << ")" << endl;
            }

            bool allOK = true;

            fvScalarMatrix ePotentialGasEqn
            (
                fvm::laplacian(epsilon_, ePotential_)
            );

            if (hasPlasma())
            {
                ePotentialGasEqn == -species().chargeDensity();
            }

            scalar resGas = ePotentialGasEqn.solve().initialResidual();
            if (resGas >= nonCoupledTolerance_) allOK = false;

            for (dielectricRegion& reg : dielectrics_)
            {
                fvScalarMatrix ePotentialDielectricEqn
                (
                    fvm::laplacian(reg.epsilon(), reg.ePotential())
                );

                scalar resDiel = dielEqn.solve().initialResidual();
                if (resDiel >= nonCoupledTolerance_) allOK = false;
            }

            if (allOK)
            {
                if (maxNonCoupledIterations_ > 1)
                {
                    Info<< ">>> ePotential converged in " << iter 
                        << " iterations." << endl;
                }
                break;
            }
            else if (iter == maxNonCoupledIterations_)
            {
                Info<< ">>> WARNING: ePotential did NOT converge after " 
                    << iter << " iterations." << endl;
            }
        }
    }

    E_ = -fvc::grad(ePotential_);
    Emag_ = mag(E_);
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();
    // for (dielectricRegion& reg : dielectrics_) reg.updateE();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
