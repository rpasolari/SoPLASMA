/*---------------------------------------------------------------------------*\
  File: threeGroupEddington.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::threeGroupEddington.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"
#include "mixedFvPatchFields.H"

#include "SoPLASMAConstants.H"
#include "threeGroupEddington.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(threeGroupEddington, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    threeGroupEddington,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

threeGroupEddington::threeGroupEddington
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
    PsiStar_j_(),
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
        boundaryConditionType_ != "marshak"
     && boundaryConditionType_ != "zero"
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Unknown boundaryConditionType '" << boundaryConditionType_
            << "'." << nl
            << "Valid options are: (marshak | zero)" << nl
            << exit(FatalIOError);
    }

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
    const List<Tuple2<scalar, scalar>> fitting
    (
        coeffs.lookup("fittingParameters")
    );

    if (fitting.size() != nGroups_)
    {
        FatalIOErrorInFunction(coeffs)
            << "The three-group Eddington model requires exactly "
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
        boundaryConditionType_ == "marshak"
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

    // Solved photon distribution components Psi*_j, one per group
    PsiStar_j_.resize(nGroups_);

    forAll(PsiStar_j_, j)
    {
        PsiStar_j_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "PsiStar_" + Foam::name(j + 1),
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

void threeGroupEddington::correct()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (!mesh_.changing() && (timeIndex % solveInterval_ != 0))
    {
        return;
    }

    Info<< "Solving photoionization (three-group Eddington, "
        << boundaryConditionType_ << " BC)..." << endl;

    const volScalarField& I =
        mesh_.lookupObject<volScalarField>(IName_);

    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    // Solve each Eddington group equation independently
    forAll(A_, j)
    {
        // 3*(lambda_j*p)^2  with units [1/m^2]
        const dimensionedScalar mu2
        (
            "mu2",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            3.0*sqr(lambda_[j] * p_.value())
        );

        // 3*lambda_j*p/c  with units [1/m]
        const dimensionedScalar kappa
        (
            "kappa",
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            3.0*lambda_[j] * p_.value() / c_.value()
        );

        const dimensionedScalar Ajpc
        (
            "Ajpc",
            dimensionSet(0, 0, -1, 0, 0, 0, 0),
            A_[j]*p_.value()*c_.value()
        );

        // Marshak BC: the mixed coefficients (valueFraction, refValue = 0,
        // refGrad = 0) depend on lambda_j, p and the mesh geometry, not on
        // the field's own value. p_ is currently a constant set at
        // construction, so it is sufficient to guard this update on
        // mesh_.changing() alone.
        // NOTE: if p_ is later promoted to a non-uniform/time-varying
        // field, this condition must also trigger whenever p_ changes, not
        // only when the mesh changes.
        //
        // Zero BC: no per-step update needed; the fixedValue patches
        // constructed with zero in the constructor remain correct
        if (boundaryConditionType_ == "marshak")
        {
            forAll(bmesh, patchi)
            {
                if (!polyPatch::constraintType(bmesh[patchi].type()))
                {
                    mixedFvPatchScalarField& pf =
                        refCast<mixedFvPatchScalarField>
                        (
                            PsiStar_j_[j].boundaryFieldRef()[patchi]
                        );

                    const scalar k = 1.5*lambda_[j]*p_.value();

                    const scalarField& deltaCoeffs =
                        pf.patch().deltaCoeffs();

                    pf.valueFraction() = k/(k + deltaCoeffs);
                }
            }
        }

        fvScalarMatrix PsiStarjEqn
        (
            fvm::laplacian(PsiStar_j_[j])
          - fvm::Sp(mu2, PsiStar_j_[j])
         ==
          - kappa*I
        );

        PsiStarjEqn.solve(mesh_.solver(PsiStar_j_[j].name()));

        // Recover the photoionization source contribution for this group:
        //   S_j = A_j * p * c * Psi*_j
        Sph_j_[j] = Ajpc*PsiStar_j_[j];
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
