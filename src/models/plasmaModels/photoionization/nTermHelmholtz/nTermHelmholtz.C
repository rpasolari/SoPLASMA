/*---------------------------------------------------------------------------*\
  File: N2O2HelmholtzPhotoionizationModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::nTermHelmholtz.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"

#include "nTermHelmholtz.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(nTermHelmholtz, 0);
addToRunTimeSelectionTable
(
    photoionizationModel,
    nTermHelmholtz,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nTermHelmholtz::nTermHelmholtz
(
    const fvMesh& mesh
)
:
    photoionizationModel(mesh),
    coefficientUnits_(word::null),
    nTerms_(0),
    p_("p", dimensionSet(1, 0, -2, 0, 0, 0, 0), 0.0),
    IName_(word::null),
    A_(),
    lambda_(),
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
        AToSI      = 1.0 / (sqr(cmToM) * sqr(torrToPa));
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

    // Fitting parameters (lambda  A) pairs, with declared count check
    const List<Tuple2<scalar, scalar>> fitting
    (
        coeffs.lookup("fittingParameters")
    );

    if (fitting.empty())
    {
        FatalIOErrorInFunction(coeffs)
            << "Empty 'fittingParameters' list." << nl
            << "At least one (lambda  A) pair is required." << nl
            << exit(FatalIOError);
    }

    nTerms_ = coeffs.get<label>("nTerms");

    if (nTerms_ != fitting.size())
    {
        FatalIOErrorInFunction(coeffs)
            << "Inconsistent term count: nTerms = " << nTerms_
            << " but fittingParameters list has " << fitting.size()
            << " entries." << nl
            << exit(FatalIOError);
    }

    A_.setSize(nTerms_);
    lambda_.setSize(nTerms_);

    forAll(fitting, j)
    {
        lambda_[j] = fitting[j].first()  * lambdaToSI;
        A_[j]      = fitting[j].second() * AToSI;
    }

    // Partial source fields, one per Helmholtz term.
    const polyBoundaryMesh& bmesh = mesh_.boundaryMesh();

    wordList patchTypes
    (
        bmesh.size(),
        fixedValueFvPatchScalarField::typeName
    );

    forAll(bmesh, patchi)
    {
        if (polyPatch::constraintType(bmesh[patchi].type()))
        {
            patchTypes[patchi] = bmesh[patchi].type();
        }
    }

    Sph_j_.resize(nTerms_);

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
                dimensionedScalar("zero", Sph_.dimensions(), 0.0),
                patchTypes
            )
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nTermHelmholtz::correct()
{
    const label timeIndex = mesh_.time().timeIndex();

    if (!mesh_.changing() && (timeIndex % solveInterval_ != 0))
    {
        return;
    }

    Info<< "Solving photoionization (" << nTerms_
        << "-term Helmholtz)..." << endl;

    const volScalarField& I =
        mesh_.lookupObject<volScalarField>(IName_);

    // Solve each Helmholtz equation independently
    forAll(A_, j)
    {
        // (lambda_j * p)^2  with units [1/m^2]
        const dimensionedScalar mu2
        (
            "mu2",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            sqr(lambda_[j] * p_.value())
        );

        // A_j * p^2  with units [1/m^2]
        const dimensionedScalar Ap2
        (
            "Ap2",
            dimensionSet(0, -2, 0, 0, 0, 0, 0),
            A_[j] * sqr(p_.value())
        );

        fvScalarMatrix SphjEqn
        (
            fvm::laplacian(Sph_j_[j])
          - fvm::Sp(mu2, Sph_j_[j])
         ==
          - Ap2*I
        );

        SphjEqn.solve(mesh_.solver(Sph_j_[j].name()));
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
