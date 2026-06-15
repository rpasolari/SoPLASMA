/*---------------------------------------------------------------------------*\
  File: singleRegionPoissonModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::singleRegionPoissonModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "singleRegionPoissonModel.H"
#include "plasmaSimulationProfiler.H"

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

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void singleRegionPoissonModel::updateDerivedFields()
{
    plasmaSimulationProfiler::start("Electromagnetics", "Calc phiE");
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc phiE");

    plasmaSimulationProfiler::start("Electromagnetics", "Calc E");
    if (EScheme_ == "reconstruct")
    {
        E_ = fvc::reconstruct(phiE_);
    }
    else // "grad"
    {
        E_ = -fvc::grad(ePotential_);
    }
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc E");

    plasmaSimulationProfiler::start("Electromagnetics", "Correct E boundary conditions");
    E_.correctBoundaryConditions();
    plasmaSimulationProfiler::stop("Electromagnetics", "Correct E boundary conditions");


    plasmaSimulationProfiler::start("Electromagnetics", "Calc Emag");
    Emag_ = mag(E_);
    Emag_.correctBoundaryConditions();
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc Emag");


    plasmaSimulationProfiler::start("Electromagnetics", "Calc reducedE");
    if (backgroundDensityFieldPtr_)
    {
        reducedE_ = Emag_ / *backgroundDensityFieldPtr_;
    }
    else
    {
        reducedE_ = Emag_ / backgroundDensityUniform_;
    }
    reducedE_.correctBoundaryConditions();
    plasmaSimulationProfiler::stop("Electromagnetics", "Calc reducedE");

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleRegionPoissonModel::singleRegionPoissonModel
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
    backgroundDensityFieldPtr_(nullptr)
{
    if (dielectricMeshes.size() > 0)
    {
        FatalErrorInFunction
            << "Dielectric regions were detected (" << dielectricMeshes.size()
            << "), but you are using the 'singleRegionPoisson' model." << nl
            << "This model ignores dielectrics and solves only in the gas."
            << nl << "Switch to 'multiRegionPoisson' model in your "
            << "properties file or remove the dielectric regions." << nl
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

    epsilonR_ = coeffs.getOrDefault<scalar>("dielectricConstant", 1.0);
    epsilon_ = dimensionedScalar
                            ("epsilon", epsilonR_ * constant::plasma::epsilon0);

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
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

//- Explicit Poisson branch
void singleRegionPoissonModel::solve()
{
    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        plasmaSimulationProfiler::start("Electromagnetics", "Build ePotentialEqn");
        fvScalarMatrix ePotentialEqn
        (
            fvm::laplacian(epsilon_, ePotential_)
         == -chargeDensity_
        );
        plasmaSimulationProfiler::stop("Electromagnetics", "Build ePotentialEqn");

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }
        plasmaSimulationProfiler::start("Electromagnetics", "Solve ePotentialEqn");
        ePotentialEqn.solve();
        plasmaSimulationProfiler::stop("Electromagnetics", "Solve ePotentialEqn");
    }

    // plasmaSimulationProfiler::start("emupdateDerivedFields");
    updateDerivedFields();
    // plasmaSimulationProfiler::stop("emupdateDerivedFields");
}

//- Semi-implicit Poisson branch
void singleRegionPoissonModel::solve
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

    const dimensionedScalar deltaT = mesh_.time().deltaT();

    const volScalarField effEps(epsilon_ + deltaT * electricalConductivity);
    const volScalarField rhsSource(-chargeDensity_ - deltaT * diffusiveChargeSource);

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix ePotentialEqn
        (
            fvm::laplacian
            (
                effEps, ePotential_
            )
         == rhsSource
        );

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }

        ePotentialEqn.solve();
    }

    updateDerivedFields();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
