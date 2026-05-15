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
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void singleRegionPoissonModel::solve()
{
    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix ePotentialEqn(PoissonEquation());

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }

        ePotentialEqn.solve();
    }

    updateDerivedFields();
}


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

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        fvScalarMatrix ePotentialEqn
        (
            PoissonEquation(electricalConductivity, diffusiveChargeSource)
        );

        if (nonOrth < nNonOrthCorr_)
        {
            ePotentialEqn.relax();
        }

        ePotentialEqn.solve();
    }

    updateDerivedFields();
}


tmp<fvScalarMatrix> singleRegionPoissonModel::PoissonLHSMatrix() const
{
    return fvm::laplacian(epsilon_, ePotential_);
}


tmp<fvScalarMatrix> singleRegionPoissonModel::PoissonEquation() const
{
    tmp<fvScalarMatrix> tEqn = PoissonLHSMatrix();

    tEqn.ref() == -chargeDensity_;

    return tEqn;
}


tmp<fvScalarMatrix> singleRegionPoissonModel::PoissonEquation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
