/*---------------------------------------------------------------------------*\
  File: ddWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ddWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "ddWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ddWallFluxImplicitFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

//- Thermal velocity for a single constant temperature (dimensionedScalar)
dimensionedScalar ddWallFluxImplicitFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T
) const
{
    return 0.25 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann * T)
        /
        (constant::mathematical::pi * m)
    );
}

//- Thermal velocity for a field temperature (scalarField)
tmp<scalarField> ddWallFluxImplicitFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const scalarField& T
) const
{
    return 0.25 * sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann.value() * T)
        /
        (constant::mathematical::pi * m.value())
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("none"),
    TValue_("T", dimTemperature, 300.0)
{
    this->refValue()      = 0.0;
    this->refGrad()       = 0.0;
    this->valueFraction() = 0.0;
}

// Dictionary Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("none"),
    TValue_("T", dimTemperature, 300.0)
{
    const entry& e = dict.lookupEntry("T", keyType::LITERAL);
    ITstream& is = e.stream();

    if (is.peek().isWord())
    {
        is >> TName_;
    }
    else if (is.peek().isNumber())
    {
        is >> TValue_.value();
        TName_ = "none";
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Entry 'T' must be either a word (field name) or a scalar value."
            << exit(FatalIOError);
    }

    if (this->readMixedEntries(dict))
    {
        // Full restart or values provided in dictionary
    }
    else
    {
        this->refValue()      = 0.0;
        this->refGrad()       = 0.0;
        this->valueFraction() = 0.0;
    }
}

// Mapping Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from another patch field)
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from patch field and new internal field)
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ddWallFluxImplicitFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const fvPatch& p = patch();
    if (p.size() == 0)
    {
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    // Access the mesh from the patch
    const fvMesh& mesh = p.boundaryMesh().mesh();
    // Access the normal vector and delta coeffs
    const vectorField nf = p.nf();
    const scalarField& delta = p.deltaCoeffs()

    // 1. Determine species name (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

    // 2. Lookup Temperature (T)
    tmp<scalarField> tT;
    if (TName_ == "none")
    {
        tT = tmp<scalarField>::New(p.size(), TValue_.value());
    }
    else if (db().foundObject<volScalarField>(TName_))
    {
        tT = p.lookupPatchField<volScalarField, scalar>(TName_);
    }
    else if (db().foundObject<UniformDimensionedField<scalar>>(TName_))
    {
        const auto& udf = db().lookupObject<UniformDimensionedField<scalar>>(TName_);
        tT = tmp<scalarField>::New(p.size(), udf.value());
    }
    else
    {
        FatalErrorInFunction << "Temperature '" << TName_ << "' not found." << exit(FatalError);
    }
    const scalarField& T = tT();

    // 3. Lookup Transport Registry and Species Data
    const plasmaTransport& transport = db().lookupObject<plasmaTransport>("plasmaTransport");
    const plasmaSpecies& speciesDB = transport.species();
    const label speciesID = speciesDB.speciesID(speciesName);
    const dimensionedScalar& m = speciesDB.speciesMass(speciesID);
    const scalar& Z = speciesDB.speciesChargeNumber(speciesID);

    // 4. Access the Model and Extract mu/D
    const plasmaTransportModel& baseModel = transport.model(speciesID);
    
    // We strictly require driftDiffusion for this BC
    if (!isA<driftDiffusion>(baseModel))
    {
        FatalErrorInFunction << "Species " << speciesName 
            << " must use driftDiffusion transport model." << exit(FatalError);
    }
    const driftDiffusion& ddModel = refCast<const driftDiffusion>(baseModel);
tmp<fvPatchScalarField> tmuPatch = 
    fvPatchField<scalar>::New("calculated", p, this->internalField());

tmp<fvPatchScalarField> tDPatch = 
    fvPatchField<scalar>::New("calculated", p, this->internalField());
    ddModel.mobility().muPatch(tmuPatch.ref(), p.index());
    ddModel.diffusivity().DPatch(tDPatch.ref(), p.index());
    const scalarField& muf = tmuPatch();
    const scalarField& Df = tDPatch();
    const fvPatchField<vector>& Ef = p.lookupPatchField<volVectorField, vector>("E");

    // 5. Physics Calculations
    word scheme = ddModel.fluxScheme();
    scalarField uDrift_n = (muf * Ef) & nf;

    // Get wall velocity from the child class (e.g., thermal flux or recombination)
    tmp<scalarField> tUWall = this->calcWallVelocity(m, T, uDrift_n, Z);
    const scalarField& uWall = tUWall();

    // 6. Set Mixed BC parameters
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->operator==(this->patchInternalField());
    scalarField& f = this->valueFraction();

    scalarField D_delta = Df * delta;

    if (scheme == "ScharfetterGummel")
    {
        auto Bern = [](scalar x) -> scalar
        {
            const scalar ax = mag(x);
            if (ax < 1e-4) return 1.0 - 0.5*x + (x*x)/12.0;
            if (x > 100.0)  return 0.0;
            if (x < -100.0) return -x;
            return x / (Foam::exp(x) - 1.0);
        };

        forAll(p, faceI)
        {
            scalar Pe = uDrift_n[faceI] / (D_delta[faceI] + VSMALL);
            scalar num = uWall[faceI] - D_delta[faceI] * Pe;
            scalar den = D_delta[faceI] * Bern(Pe) + uWall[faceI];
            f[faceI] = num / (den + VSMALL);
        }
    }
    else // Standard Upwind/Central scheme
    {
        forAll(p, faceI)
        {
            f[faceI] = uWall[faceI] / (D_delta[faceI] + uWall[faceI] + SMALL);
        }
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}

void ddWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    // 1. Write standard Mixed BC entries (valueFraction, refValue, etc.)
    mixedFvPatchScalarField::write(os);

    // 2. Write our custom entries so the simulation can be restarted
    if (TName_ != "none")
    {
        os.writeEntry("T", TName_);
    }
    else
    {
        os.writeEntry("T", TValue_.value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
