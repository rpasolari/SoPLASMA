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

singleRegionPoissonModel::singleRegionPoissonModel
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
    transportPtr_(nullptr)
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

    if (PoissonScheme_ != "explicit" && PoissonScheme_ != "semiImplicit")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown PoissonScheme '" << PoissonScheme_ << "'." << nl
            << "Valid options are: (explicit | semiImplicit)" << nl
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void singleRegionPoissonModel::solve()
{
    if (PoissonScheme_ == "semiImplicit" && transportPtr_)
    {
        this->solve(*transportPtr_);
    }
    else
    {
        fvScalarMatrix ePotentialEqn
        (
            fvm::laplacian(epsilon_, ePotential_)
        );

        if (hasPlasma())
        {
            ePotentialEqn == -species().chargeDensity();
        }

        ePotentialEqn.solve();

        E_ = -fvc::grad(ePotential_);
        Emag_ = mag(E_);
        phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();
    }
}

void singleRegionPoissonModel::solve(const plasmaTransportModel& transport)
{
    transportPtr_ = &transport;

    if (!hasPlasma())
    {
        this->solve();
        return;
    }

    const scalar deltaT = mesh_.time().deltaTValue();
    const volScalarField& chargeDensity = species().chargeDensity();

    fvScalarMatrix ePotentialEqn
    (
        fvm::laplacian
        (
            epsilon_ + deltaT * transport.electricalConductivity(),
            ePotential_
        )
     ==
        -chargeDensity - deltaT * transport.diffusiveChargeSource()
    );

    ePotentialEqn.solve();

    E_ = -fvc::grad(ePotential_);
    Emag_ = mag(E_);
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();
}

tmp<fvScalarMatrix> singleRegionPoissonModel::PoissonMatrix() const
{
    return fvm::laplacian(epsilon_, ePotential_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
