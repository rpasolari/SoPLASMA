/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dielectricRegion
    
Description
    A container class for dielectric region fields and properties.

    It is designed to be used in multi-region Poisson solvers to eliminate the
    need of multiple synchronized pointer lists.

    Holds per-region storage for ePotential, E, phiE, and permittivity. Each
    region uses a single uniform permittivity (no spatial variation), set from
    'dielectricConstant' in the region's sub-dictionary.

\*---------------------------------------------------------------------------*/

#ifndef dielectricRegion_H
#define dielectricRegion_H

#include "fvCFD.H"

#include "plasmaConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dielectricRegion Declaration
\*---------------------------------------------------------------------------*/

class dielectricRegion
{

protected:
   
    //- Reference to the mesh for this specific dielectric region
    const fvMesh& regionMesh_;

    //- The dictionary containing the information for the specific region
    dictionary regionDict_;

    //- Electric potential field [V]
    volScalarField ePotential_;

    //- Electric field [V/m]
    volVectorField E_;

    //- Flux of the electric field [V * m2]
    surfaceScalarField phiE_;

    //- Surface charge density (owned by dielectric) [C/m2]
    volScalarField surfCharge_;

    //- Dielectric permittivity [F/m]
    dimensionedScalar epsilon_;

    //- Relative permittivity of the medium
    scalar epsilonR_;

    //- Numerical scheme for computing E
    word EScheme_;

public:

    // Constructors

        //- Construct from mesh and dictionary
        dielectricRegion(const fvMesh& regionMesh, const dictionary& regionDict)
        :
            regionMesh_(regionMesh),
            regionDict_(regionDict),
            ePotential_
            (
                IOobject
                (
                    "ePotential",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_
            ),
            E_
            (
                IOobject
                (
                    "E",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(1, 1, -3, 0, 0, -1, 0),
                    vector::zero
                )
            ),
            phiE_
            (
                IOobject
                (
                    "phiE",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -fvc::snGrad(ePotential_) * regionMesh_.magSf()
            ),
            surfCharge_
            (
                IOobject
                (
                    "surfCharge",
                    regionMesh_.time().timeName(),
                    regionMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                regionMesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, -2, 1, 0, 0, 1, 0),
                    0.0
                )
            ),
            epsilon_
            (
                "epsilon",
                dimensionSet(-1, -3, 4, 0, 0, 2, 0),
                0.0
            ),
            epsilonR_
            (
                regionDict_.getOrDefault<scalar>("dielectricConstant", 1.0)
            ),
            EScheme_
            (
                regionDict_.getOrDefault<word>("EScheme", "reconstruct")
            )
        {
            epsilon_ = epsilonR_ * constant::plasma::epsilon0;

            if (EScheme_ != "grad" && EScheme_ != "reconstruct")
            {
                FatalIOErrorInFunction(regionDict_)
                    << "Unknown EScheme '" << EScheme_ << "' in region '"
                    << regionMesh_.name() << "'." << nl
                    << "Valid options are: (grad | reconstruct)" << nl
                    << exit(FatalIOError);
            }

            // Register a lightweight electricProperties dictionary in this
            // region's mesh registry
            IOdictionary* epsilonDict = new IOdictionary
            (
                IOobject
                (
                    "electricProperties",
                    regionMesh_.time().constant(),
                    regionMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            );
            epsilonDict->add("dielectricConstant", epsilonR_);
            regionMesh_.thisDb().store(epsilonDict);
        }

    //- Destructor
    virtual ~dielectricRegion() = default;

    // Member functions

        //- Update the electric field
        void updateE()
        {
            phiE_ = -fvc::snGrad(ePotential_) * regionMesh_.magSf();

            if (EScheme_ == "reconstruct")
            {
                E_ = fvc::reconstruct(phiE_);
            }
            else // "grad"
            {
                E_ = -fvc::grad(ePotential_);
            }
            E_.correctBoundaryConditions();
        }

        //- Return the LHS of the Poisson equation
        tmp<fvScalarMatrix> PoissonMatrix() const
        {
            return fvm::laplacian(epsilon_, ePotential_);
        }

        tmp<fvScalarMatrix> PoissonEquation() const
        {
            // In a dielectric, the Poisson equation reduces to Laplace
            // (no free charges in the bulk; surface charges enter via BCs)
            return PoissonMatrix();
        }

    // Access

        inline const fvMesh& mesh() const { return regionMesh_; }
        inline const dictionary& dict() const { return regionDict_; }
        inline const volScalarField& ePotential() const { return ePotential_; }
        inline volScalarField& ePotential() { return ePotential_; }
        inline const volVectorField& E() const { return E_; }
        inline const surfaceScalarField& phiE() const { return phiE_; }
        inline const volScalarField& surfCharge() const { return surfCharge_; }
        inline volScalarField& surfCharge() { return surfCharge_; }
        inline const dimensionedScalar& epsilon() const { return epsilon_; }
        inline scalar epsilonR() const { return epsilonR_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

