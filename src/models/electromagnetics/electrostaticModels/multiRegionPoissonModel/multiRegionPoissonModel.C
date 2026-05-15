/*---------------------------------------------------------------------------*\
  File: multiRegionPoissonModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::multiRegionPoissonModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
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

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void multiRegionPoissonModel::updateDerivedFields()
{
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();

    if (EScheme_ == "reconstruct")
    {
        E_ = fvc::reconstruct(phiE_);
    }
    else // "grad"
    {
        E_ = -fvc::grad(ePotential_);
    }
    E_.correctBoundaryConditions();

    Emag_ = mag(E_);
    Emag_.correctBoundaryConditions();

    if (backgroundDensityFieldPtr_)
    {
        reducedE_ = Emag_ /
            (
                *backgroundDensityFieldPtr_
              + dimensionedScalar
                (
                    "s",
                    backgroundDensityUniform_.dimensions(),
                    SMALL
                )
            );
    }
    else
    {
        reducedE_ = Emag_ /
            (
                backgroundDensityUniform_
              + dimensionedScalar
                (
                    "s",
                    backgroundDensityUniform_.dimensions(),
                    SMALL
                )
            );
    }
    reducedE_.correctBoundaryConditions();

    for (dielectricRegion& reg : dielectrics_)
    {
        reg.updateE();
    }
}

void multiRegionPoissonModel::solveCoupled()
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(monolithically)" << endl;

    ePotential_.correctBoundaryConditions();
    for (dielectricRegion& reg : dielectrics_)
    {
        reg.ePotential().correctBoundaryConditions();
    }

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix gasEqn(PoissonEquation());

        fvMatrixAssemblyPtr_->addFvMatrix(gasEqn);

        for (dielectricRegion& reg : dielectrics_)
        {
            fvMatrixAssemblyPtr_->addFvMatrix(reg.PoissonMatrix().ref());
        }

        fvMatrixAssemblyPtr_->solve();

        ePotential_.correctBoundaryConditions();
        for (dielectricRegion& reg : dielectrics_)
        {
            reg.ePotential().correctBoundaryConditions();
        }

        fvMatrixAssemblyPtr_->clear();
    }
}

void multiRegionPoissonModel::solveCoupled
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(monolithically)" << endl;

    ePotential_.correctBoundaryConditions();
    for (dielectricRegion& reg : dielectrics_)
    {
        reg.ePotential().correctBoundaryConditions();
    }

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix gasEqn
        (
            PoissonEquation(electricalConductivity, diffusiveChargeSource)
        );

        fvMatrixAssemblyPtr_->addFvMatrix(gasEqn);

        for (dielectricRegion& reg : dielectrics_)
        {
            fvMatrixAssemblyPtr_->addFvMatrix(reg.PoissonMatrix().ref());
        }

        fvMatrixAssemblyPtr_->solve();

        ePotential_.correctBoundaryConditions();
        for (dielectricRegion& reg : dielectrics_)
        {
            reg.ePotential().correctBoundaryConditions();
        }

        fvMatrixAssemblyPtr_->clear();
    }
}

void multiRegionPoissonModel::solveSegregated()
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(segregated)" << endl;

    for (int iter = 1; iter <= maxNonCoupledIterations_; ++iter)
    {
        if (maxNonCoupledIterations_ > 1)
        {
            Info<< "  -ePotential iteration (" << iter << "/"
                << maxNonCoupledIterations_ << ")" << endl;
        }

        bool allOK = true;

        {
            scalar resGas = 0.0;

            for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
            {
                fvScalarMatrix gasEqn(PoissonEquation());

                if (nonOrth < nNonOrthCorr_)
                    gasEqn.relax();

                Info<< "    -Solving for ePotential (gas region: "
                    << mesh_.name() << ")" << endl;

                auto sp = gasEqn.solve();

                if (nonOrth == 0)
                    resGas = sp.initialResidual();

                ePotential_.correctBoundaryConditions();
            }

            if (resGas >= nonCoupledTolerance_) allOK = false;
        }

        for (dielectricRegion& reg : dielectrics_)
        {
            scalar resDiel = 0.0;

            reg.ePotential().correctBoundaryConditions();

            for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
            {
                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(reg.epsilon(), reg.ePotential())
                );

                if (nonOrth < nNonOrthCorr_)
                    dielEqn.relax();

                Info<< "    -Solving for ePotential (dielectric: "
                    << reg.mesh().name() << ")" << endl;

                auto sp = dielEqn.solve();

                if (nonOrth == 0)
                    resDiel = sp.initialResidual();

                reg.ePotential().correctBoundaryConditions();
            }

            if (resDiel >= nonCoupledTolerance_) allOK = false;
        }

        ePotential_.correctBoundaryConditions();

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

void multiRegionPoissonModel::solveSegregated
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    Info<< "Solving for ePotential in non-coupled regions "
        << "(segregated)" << endl;

    for (int iter = 1; iter <= maxNonCoupledIterations_; ++iter)
    {
        if (maxNonCoupledIterations_ > 1)
        {
            Info<< "  -ePotential iteration (" << iter << "/"
                << maxNonCoupledIterations_ << ")" << endl;
        }

        bool allOK = true;

        for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
        {
            // Rebuild augmented gas equation each pass
            fvScalarMatrix gasEqn
            (
                PoissonEquation(electricalConductivity, diffusiveChargeSource)
            );

            if (nonOrth < nNonOrthCorr_)
            {
                gasEqn.relax();
            }

            Info<< "    -Solving for ePotential (gas region: "
                << mesh_.name() << ")" << endl;

            scalar resGas = gasEqn.solve().initialResidual();

            ePotential_.correctBoundaryConditions();

            if (nonOrth == nNonOrthCorr_)
            {
                if (resGas >= nonCoupledTolerance_) allOK = false;
            }

            for (dielectricRegion& reg : dielectrics_)
            {
                // Dielectric always uses Laplace (no source term)
                fvScalarMatrix dielEqn
                (
                    fvm::laplacian(reg.epsilon(), reg.ePotential())
                );

                if (nonOrth < nNonOrthCorr_)
                {
                    dielEqn.relax();
                }

                Info<< "    -Solving for ePotential (dielectric: "
                    << reg.mesh().name() << ")" << endl;

                scalar resDiel = dielEqn.solve().initialResidual();

                reg.ePotential().correctBoundaryConditions();

                if (nonOrth == nNonOrthCorr_)
                {
                    if (resDiel >= nonCoupledTolerance_) allOK = false;
                }
            }
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
multiRegionPoissonModel::multiRegionPoissonModel
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
:
    electromagneticsModel(mesh, dielectricMeshes),
    EScheme_("reconstruct"),
    PoissonScheme_("explicit"),
    nNonOrthCorr_(0),
    backgroundDensityUniform_
    (
        "backgroundDensity",
        dimensionSet(0, -3, 0, 0, 0, 0, 0),
        0.0
    ),
    backgroundDensityFieldPtr_(nullptr),
    dielectrics_(dielectricMeshes.size()),
    coupled_(false),
    fvMatrixAssemblyPtr_(nullptr),
    maxNonCoupledIterations_(100),
    nonCoupledTolerance_(1e-6)
{
    if (dielectricMeshes.size() == 0)
    {
        FatalErrorInFunction
            << "No dielectric regions were detected, but you are using the "
            << "'multiRegionPoisson' model." << nl
            << "This model requires at least one dielectric mesh." << nl
            << "Switch to 'singleRegionPoisson' model in your properties "
            << "file if you only have a gas region." << nl
            << exit(FatalError);
    }

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

    EScheme_ = coeffs.getOrDefault<word>("EScheme", "reconstruct");
    if (EScheme_ != "grad" && EScheme_ != "reconstruct")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown EScheme '" << EScheme_ << "'." << nl
            << "Valid options are: (grad | reconstruct)" << nl
            << exit(FatalIOError);
    }

    PoissonScheme_ = coeffs.getOrDefault<word>("PoissonScheme", "explicit");
    if (PoissonScheme_ != "explicit" && PoissonScheme_ != "semiImplicit")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown PoissonScheme '" << PoissonScheme_ << "'." << nl
            << "Valid options are: (explicit | semiImplicit)" << nl
            << exit(FatalIOError);
    }

    nNonOrthCorr_ = coeffs.getOrDefault<label>("nNonOrthogonalCorrectors", 0);

    backgroundDensityUniform_.value() =
        coeffs.getOrDefault<scalar>("backgroundDensity", 2.5e25);

    const dictionary& nonCoupledResidualControl =
        coeffs.subOrEmptyDict("nonCoupledResidualControl");

    maxNonCoupledIterations_ =
        nonCoupledResidualControl.getOrDefault<label>("maxIter", 100);

    nonCoupledTolerance_ =
        nonCoupledResidualControl.getOrDefault<scalar>("tolerance", 1e-6);

    // Build dielectric regions
    forAll(dielectricMeshes_, i)
    {
        const word& regionName = dielectricMeshes_[i].name();
        if (!coeffs.found(regionName))
        {
            FatalIOErrorInFunction(coeffs)
                << "Region '" << regionName << "' not found in "
                << coeffsName << " dictionary." << nl
                << "Every dielectric mesh must have a corresponding "
                << "sub-dictionary." << exit(FatalIOError);
        }
        const dictionary& regionDict = coeffs.subDict(regionName);

        dielectrics_.set
        (
            i,
            new dielectricRegion(dielectricMeshes_[i], regionDict)
        );
    }

    // Detect implicit coupling
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
        Info<< "    Implicit coupling detected." << endl;

        // Validate all mapped patches are consistent
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
            << "SEGREGATED way" << nl << endl;
    }



    // --- DEBUG REGISTRY INSPECTION ---
Info<< nl << "/*-----------------------------------------------------------*\\" << endl;
Info<< "  DEBUG: multiRegionPoissonModel Registry Final Check" << endl;

// 1. Check the Primary Gas Mesh Registry
Info<< "  Gas Mesh (" << mesh_.name() << ") Registry:" << nl
    << "  " << mesh_.thisDb().sortedToc() << nl << endl;

// Explicitly check for the dictionary our model represents
if (mesh_.thisDb().foundObject<IOdictionary>("electromagneticsProperties"))
{
    const IOdictionary& d = mesh_.thisDb().lookupObject<IOdictionary>("electromagneticsProperties");
    Info<< "  FOUND: 'electromagneticsProperties' in gas mesh." << nl
        << "  Address: " << &d << " | Class: " << d.type() << endl;
}
else
{
    Info<< "  CRITICAL: 'electromagneticsProperties' NOT FOUND in gas mesh registry!" << endl;
}

// 2. Check each Dielectric Mesh Registry
forAll(dielectrics_, i)
{
    const fvMesh& dMesh = dielectrics_[i].mesh();
    Info<< nl << "  Dielectric Mesh (" << dMesh.name() << ") Registry:" << nl
        << "  " << dMesh.thisDb().sortedToc() << endl;

    if (dMesh.thisDb().foundObject<IOdictionary>("electricProperties"))
    {
        Info<< "  FOUND: 'electricProperties' in " << dMesh.name() << endl;
    }
    else
    {
        Info<< "  WARNING: 'electricProperties' NOT FOUND in " << dMesh.name() << endl;
    }
}

Info<< "\\*-----------------------------------------------------------*/" << nl << endl;
// --- END DEBUG REGISTRY INSPECTION ---
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void multiRegionPoissonModel::solve()
{
    if (coupled_)
    {
        solveCoupled();
    }
    else
    {
        solveSegregated();
    }

    updateDerivedFields();
}


void multiRegionPoissonModel::solve
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
)
{
    if (PoissonScheme_ == "explicit")
    {
        this->solve();
        return;
    }

    if (coupled_)
    {
        solveCoupled(electricalConductivity, diffusiveChargeSource);
    }
    else
    {
        solveSegregated(electricalConductivity, diffusiveChargeSource);
    }

    updateDerivedFields();
}


tmp<fvScalarMatrix> multiRegionPoissonModel::PoissonLHSMatrix() const
{
    return fvm::laplacian(epsilon_, ePotential_);
}


tmp<fvScalarMatrix> multiRegionPoissonModel::PoissonEquation() const
{
    tmp<fvScalarMatrix> tEqn = PoissonLHSMatrix();

    tEqn.ref() == -chargeDensity_;

    return tEqn;
}


tmp<fvScalarMatrix> multiRegionPoissonModel::PoissonEquation
(
    const volScalarField& electricalConductivity,
    const volScalarField& diffusiveChargeSource
) const
{
    const scalar deltaT = mesh_.time().deltaTValue();

    tmp<fvScalarMatrix> tEqn = fvm::laplacian
    (
        epsilon_ + deltaT * electricalConductivity,
        ePotential_
    );

    tEqn.ref() == -chargeDensity_ - deltaT * diffusiveChargeSource;

    return tEqn;
}


const dielectricRegion& multiRegionPoissonModel::dielectric
(
    const word& name
) const
{
    forAll(dielectrics_, i)
    {
        if (dielectrics_[i].mesh().name() == name)
        {
            return dielectrics_[i];
        }
    }

    wordList availableNames(dielectrics_.size());
    forAll(dielectrics_, i)
    {
        availableNames[i] = dielectrics_[i].mesh().name();
    }

    FatalErrorInFunction
        << "Dielectric region '" << name << "' not found." << nl
        << "Available: " << availableNames << nl
        << exit(FatalError);

    return dielectrics_[0]; // unreachable
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
