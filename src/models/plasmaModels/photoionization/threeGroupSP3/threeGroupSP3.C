/*---------------------------------------------------------------------------*\
  File: threeGroupSP3.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::threeGroupSP3.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "mixedFvPatchFields.H"

#include "SoPLASMAConstants.H"
#include "threeGroupSP3.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(threeGroupSP3, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    threeGroupSP3,
    dictionary
);

// * * * * * * * * * * * * * * SP3 Closure Constants  * * * * * * * * * * * //

// kappa1^2 = 3/7 - (2/7)*sqrt(6/5)
const scalar threeGroupSP3::kappa1Sqr_ = 3.0/7.0 - (2.0/7.0)*Foam::sqrt(6.0/5.0);

// kappa2^2 = 3/7 + (2/7)*sqrt(6/5)
const scalar threeGroupSP3::kappa2Sqr_ = 3.0/7.0 + (2.0/7.0)*Foam::sqrt(6.0/5.0);

// gamma1 = (5/7)*[1 - 3*sqrt(6/5)]
const scalar threeGroupSP3::gamma1_ = (5.0/7.0)*(1.0 - 3.0*Foam::sqrt(6.0/5.0));

// gamma2 = (5/7)*[1 + 3*sqrt(6/5)]
const scalar threeGroupSP3::gamma2_ = (5.0/7.0)*(1.0 + 3.0*Foam::sqrt(6.0/5.0));

// * * * * * * * * * * * * * * Larsen BC Constants  * * * * * * * * * * * * //

// alpha1 = (5/96)*(34 + 11*sqrt(6/5))
const scalar threeGroupSP3::alpha1_ =
    (5.0/96.0)*(34.0 + 11.0*Foam::sqrt(6.0/5.0));

// alpha2 = (5/96)*(34 - 11*sqrt(6/5))
const scalar threeGroupSP3::alpha2_ =
    (5.0/96.0)*(34.0 - 11.0*Foam::sqrt(6.0/5.0));

// beta1 = (5/96)*(2 - sqrt(6/5))
const scalar threeGroupSP3::beta1_ =
    (5.0/96.0)*(2.0 - Foam::sqrt(6.0/5.0));

// beta2 = (5/96)*(2 + sqrt(6/5))
const scalar threeGroupSP3::beta2_ =
    (5.0/96.0)*(2.0 + Foam::sqrt(6.0/5.0));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

threeGroupSP3::threeGroupSP3
(
    const fvMesh& mesh
)
:
    photoionizationModel(mesh),
    coefficientUnits_(word::null),
    boundaryConditionType_(word::null),
    p_("p", dimensionSet(1, 0, -2, 0, 0, 0, 0), 0.0),
    IName_(word::null),
    A_(),
    lambda_(),
    c_
    ("c", dimensionSet(0, 1, -1, 0, 0, 0, 0), constant::plasma::cLight.value()),
    larsenMaxIters_(1),
    larsenTol_(1e-6),
    phi1_j_(),
    phi2_j_(),
    Sph_j_()
{
    const word coeffsName(type() + "Coeffs");

    if (!found(coeffsName))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary '" << coeffsName << "' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& coeffs(subDict(coeffsName));

    // Boundary condition type
    boundaryConditionType_ =
        coeffs.getOrDefault<word>("boundaryConditionType", "marshak");

    if
    (
        boundaryConditionType_ != "larsen"
     && boundaryConditionType_ != "zero"
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown boundaryConditionType '" << boundaryConditionType_
            << "'." << nl
            << "Valid options are: (larsen | zero)" << nl
            << exit(FatalIOError);
    }

    larsenMaxIters_ =
        coeffs.getOrDefault<label>("larsenMaxIters", 1);

    larsenTol_ =
        coeffs.getOrDefault<scalar>("larsenTol", 1e-6);

    // Unit convention and conversion factors to SI
    coefficientUnits_ = coeffs.getOrDefault<word>("coefficientUnits", "SI");

    scalar pToSI = 1.0;
    scalar lambdaToSI = 1.0;
    scalar AToSI = 1.0;

    if (coefficientUnits_ == "TorrCm")
    {
        const scalar torrToPa = 101325.0/760.0;
        const scalar cmToM    = 0.01;

        pToSI      = torrToPa;
        lambdaToSI = 1.0 / (cmToM * torrToPa);
        AToSI      = 1.0 / (cmToM * torrToPa);
    }
    else if (coefficientUnits_ != "SI")
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown coefficientUnits '" << coefficientUnits_ << "'." << nl
            << "Valid options are: (SI | TorrCm)" << nl
            << exit(FatalIOError);
    }

    // Gas pressure (read in user units, store internally in Pa)
    const scalar pRaw = coeffs.get<scalar>("p");

    p_.value() = pRaw * pToSI;

    // Ionization source field (lookup name, verify registration)
    IName_ = coeffs.get<word>("I");

    if (!mesh_.foundObject<volScalarField>(IName_))
    {
        FatalIOErrorInFunction(coeffs)
            << "Ionization source field '" << IName_
            << "' not found in the mesh object registry." << nl
            << "It must be constructed and registered before "
            << "photoionizationModel::New() is called." << nl
            << exit(FatalIOError);
    }

    // Fitting parameters (lambda  A) pairs, fixed to the three-group model
    // (same parametrization as the three-group Eddington model)
    const List<Tuple2<scalar, scalar>> fitting
    (
        coeffs.lookup("fittingParameters")
    );

    if (fitting.size() != nGroups_)
    {
        FatalIOErrorInFunction(coeffs)
            << "The three-group SP3 model requires exactly "
            << nGroups_ << " (lambda  A) pairs in 'fittingParameters', "
            << "but " << fitting.size() << " were given." << nl
            << exit(FatalIOError);
    }

    A_.setSize(nGroups_);
    lambda_.setSize(nGroups_);

    forAll(fitting, j)
    {
        lambda_[j] = fitting[j].first()  * lambdaToSI;
        A_[j]      = fitting[j].second() * AToSI;
    }

    // Boundary patch types: mixed (Robin) on non-constraint patches to
    // support the Marshak boundary condition; constraint types (e.g.
    // empty, symmetry, cyclic) are preserved as-is.
    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    const word nonConstraintPatchType =
    (
        boundaryConditionType_ == "larsen"
      ? mixedFvPatchScalarField::typeName
      : fixedValueFvPatchScalarField::typeName
    );

    wordList patchTypes(bmesh.size(), nonConstraintPatchType);

    forAll(bmesh, patchi)
    {
        if (polyPatch::constraintType(bmesh[patchi].type()))
        {
            patchTypes[patchi] = bmesh[patchi].type();
        }
    }

    // Auxiliary photon distribution components phi1_j and phi2_j, one
    // pair per group
    phi1_j_.resize(nGroups_);
    phi2_j_.resize(nGroups_);

    forAll(phi1_j_, j)
    {
        phi1_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "phi1_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero", dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0
                ),
                patchTypes
            )
        );

        phi2_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "phi2_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero", dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0
                ),
                patchTypes
            )
        );
    }

    // Partial source fields S_j, one per group
    Sph_j_.resize(nGroups_);

    forAll(Sph_j_, j)
    {
        Sph_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "Sph_" + Foam::name(j + 1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", Sph_.dimensions(), 0.0)
            )
        );
    }
    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void threeGroupSP3::correct()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (!mesh_.changing() && (timeIndex % solveInterval_ != 0))
    {
        return;
    }

    Info<< "Solving photoionization (three-group SP3, "
        << boundaryConditionType_ << " BC)..." << endl;

    const volScalarField& I =
        mesh_.lookupObject<volScalarField>(IName_);

    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    // Solve each pair of SP3 sub-equations independently
    forAll(A_, j)
    {
        // (lambda_j*p)^2/kappa1^2  with units [1/m^2]
        const dimensionedScalar mu1Sqr
        (
            "mu1Sqr",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            sqr(lambda_[j]*p_.value())/kappa1Sqr_
        );

        // (lambda_j*p)^2/kappa2^2  with units [1/m^2]
        const dimensionedScalar mu2Sqr
        (
            "mu2Sqr",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            sqr(lambda_[j]*p_.value())/kappa2Sqr_
        );

        // (lambda_j*p)/(kappa1^2*c)  with units [s/m^2]
        const dimensionedScalar kappaSrc1
        (
            "kappaSrc1",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            (lambda_[j]*p_.value())/(kappa1Sqr_*c_.value())
        );

        // (lambda_j*p)/(kappa2^2*c)  with units [s/m^2]
        const dimensionedScalar kappaSrc2
        (
            "kappaSrc2",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            (lambda_[j]*p_.value())/(kappa2Sqr_*c_.value())
        );

        // The Larsen mixed BC diagonal coefficients (valueFraction = k_alpha
        // / (k_alpha + deltaCoeffs), refValue = 0, refGrad = 0 for diagonal-
        // only iteration 0) depend on lambda_j, p and the mesh geometry, not
        // on the field's own value. Applied here every correct() call so the
        // constructor remains BC-agnostic.
        //
        // Zero BC: no per-step update needed; the fixedValue patches
        // constructed with zero in the constructor remain correct.
        if (boundaryConditionType_ == "larsen")
        {
            forAll(bmesh, patchi)
            {
                if (!polyPatch::constraintType(bmesh[patchi].type()))
                {
                    mixedFvPatchScalarField& pf1 =
                        refCast<mixedFvPatchScalarField>
                        (
                            phi1_j_[j].boundaryFieldRef()[patchi]
                        );

                    mixedFvPatchScalarField& pf2 =
                        refCast<mixedFvPatchScalarField>
                        (
                            phi2_j_[j].boundaryFieldRef()[patchi]
                        );

                    const scalar k1 = lambda_[j]*p_.value()*alpha1_;
                    const scalar k2 = lambda_[j]*p_.value()*alpha2_;

                    const scalarField& deltaCoeffs =
                        pf1.patch().deltaCoeffs();

                    pf1.refValue() = 0.0;
                    pf1.refGrad()  = 0.0;
                    pf1.valueFraction() = k1/(k1 + deltaCoeffs);

                    pf2.refValue() = 0.0;
                    pf2.refGrad()  = 0.0;
                    pf2.valueFraction() = k2/(k2 + deltaCoeffs);
                }
            }
        }

        // Inner Picard iteration on the Larsen off-diagonal beta
        // cross-coupling.
        List<scalarField> phi1BoundaryPrev(bmesh.size());
        List<scalarField> phi2BoundaryPrev(bmesh.size());

        for (label iter = 0; iter < larsenMaxIters_; ++iter)
        {
            if (boundaryConditionType_ == "larsen" && larsenMaxIters_ > 1)
            {
                Info<< "    SP3 group " << j + 1 << ", Larsen iter "
                    << iter + 1 << "/" << larsenMaxIters_ << endl;
            }

            // Update lagged beta cross-coupling in refGrad (only for
            // larsen BC and only when iterating beyond the diagonal pass).
            if (boundaryConditionType_ == "larsen" && iter > 0)
            {
                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        mixedFvPatchScalarField& pf1 =
                            refCast<mixedFvPatchScalarField>
                            (
                                phi1_j_[j].boundaryFieldRef()[patchi]
                            );

                        mixedFvPatchScalarField& pf2 =
                            refCast<mixedFvPatchScalarField>
                            (
                                phi2_j_[j].boundaryFieldRef()[patchi]
                            );

                        pf1.refGrad() =
                            -lambda_[j]*p_.value()*beta2_
                           * phi2_j_[j].boundaryField()[patchi];

                        pf2.refGrad() =
                            -lambda_[j]*p_.value()*beta1_
                           * phi1_j_[j].boundaryField()[patchi];
                    }
                }
            }

            // Snapshot previous boundary values for convergence check,
            // one entry per non-constraint patch.
            if (boundaryConditionType_ == "larsen" && larsenMaxIters_ > 1)
            {
                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        phi1BoundaryPrev[patchi] =
                            phi1_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                        phi2BoundaryPrev[patchi] =
                            phi2_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                    }
                }
            }

            fvScalarMatrix phi1jEqn
            (
                fvm::laplacian(phi1_j_[j])
              - fvm::Sp(mu1Sqr, phi1_j_[j])
             ==
              - kappaSrc1*I
            );

            phi1jEqn.solve(mesh_.solver(phi1_j_[j].name()));

            fvScalarMatrix phi2jEqn
            (
                fvm::laplacian(phi2_j_[j])
              - fvm::Sp(mu2Sqr, phi2_j_[j])
             ==
              - kappaSrc2*I
            );

            phi2jEqn.solve(mesh_.solver(phi2_j_[j].name()));

            // Convergence check on boundary values of phi1, phi2, taking the 
            // max relative change across all non-constraint patches. Normalized
            // against the GLOBAL max of the field.
            if (boundaryConditionType_ == "larsen" && iter > 0)
            {
                const scalar globalRefPhi1 = 
                            max(mag(phi1_j_[j].primitiveField())) + SMALL;
                const scalar globalRefPhi2 = 
                            max(mag(phi2_j_[j].primitiveField())) + SMALL;

                scalar relDelta1 = 0.0;
                scalar relDelta2 = 0.0;

                forAll(bmesh, patchi)
                {
                    if (!polyPatch::constraintType(bmesh[patchi].type()))
                    {
                        const scalarField phi1New =
                            phi1_j_[j].boundaryField()[patchi]
                           .patchInternalField();
                        const scalarField phi2New =
                            phi2_j_[j].boundaryField()[patchi]
                           .patchInternalField();

                        relDelta1 = max
                        (
                            relDelta1,
                            max(mag(phi1New - phi1BoundaryPrev[patchi]))
                           /globalRefPhi1
                        );
                        relDelta2 = max
                        (
                            relDelta2,
                            max(mag(phi2New - phi2BoundaryPrev[patchi]))
                           /globalRefPhi2
                        );
                    }
                }
            }
        }

        // A_j * p * c, with units [1/s]
        const dimensionedScalar Ajpc
        (
            "Ajpc",
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
            A_[j]*p_.value()*c_.value()
        );

        // SP3 closure: Psi_j = (gamma2*phi1_j - gamma1*phi2_j)/(gamma2-gamma1)
        // Recover the photoionization source contribution for this group:
        //   S_j = A_j * p * c * Psi_j
        Sph_j_[j] =
            Ajpc
           *(gamma2_*phi1_j_[j] - gamma1_*phi2_j_[j])/(gamma2_ - gamma1_);
    }

    // Sum partial sources into the total photoionization source field
    Sph_ == dimensionedScalar(Sph_.dimensions(), Zero);

    forAll(Sph_j_, j)
    {
        Sph_ += Sph_j_[j];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
