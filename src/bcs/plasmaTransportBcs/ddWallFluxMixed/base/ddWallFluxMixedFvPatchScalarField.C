/*---------------------------------------------------------------------------*\
  File: ddWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ddWallFluxMixedFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "ddWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ddWallFluxMixedFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

//- Thermal velocity for a single constant temperature (dimensionedScalar)
dimensionedScalar ddWallFluxMixedFvPatchScalarField::calcThermalVelocity
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
tmp<scalarField> ddWallFluxMixedFvPatchScalarField::calcThermalVelocity
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
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
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
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
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
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf,
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
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from patch field and new internal field)
ddWallFluxMixedFvPatchScalarField::ddWallFluxMixedFvPatchScalarField
(
    const ddWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ddWallFluxMixedFvPatchScalarField::updateCoeffs()
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
    const scalarField& delta = p.deltaCoeffs();

    // Determine species name (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

    // Lookup Temperature (T)
    tmp<scalarField> tT;
    if (TName_ == "none")
    {
        tT = tmp<scalarField>::New(p.size(), TValue_.value());
    }
    else if (db().foundObject<volScalarField>(TName_))
    {
        tT = p.lookupPatchField<volScalarField, scalar>(TName_);
    }
    else
    {
        FatalErrorInFunction
            << "Temperature field '" << TName_ << "' not found in registry." 
            << nl << "Either set TName to 'none' and provide a constant TValue,"
            << " or ensure the field exists in the mesh registry." << nl
            << exit(FatalError);
    }
    const scalarField& T = tT();

    // Lookup Transport Registry and Species Data
    if (!db().foundObject<plasmaTransport>("plasmaTransport"))
    {
        FatalErrorInFunction
            << "plasmaTransport not found in registry." << nl
            << exit(FatalError);
    }
    const plasmaTransport& transport =
                          db().lookupObject<plasmaTransport>("plasmaTransport");

    const plasmaSpecies& speciesDB = transport.species();

    const label speciesID = speciesDB.speciesID(speciesName);
    const dimensionedScalar& m = speciesDB.speciesMass(speciesID);
    const scalar Z = speciesDB.speciesChargeNumber(speciesID);

    // Access the drift-diffusion model
    const plasmaTransportModel& baseModel = transport.model(speciesID);

    if (!isA<driftDiffusion>(baseModel))
    {
        FatalErrorInFunction
            << "Species '" << speciesName << "' must use the driftDiffusion "
            << "transport model for this boundary condition." << nl
            << "Current model: " << baseModel.type() << nl
            << exit(FatalError);
    }

    const driftDiffusion& ddModel = refCast<const driftDiffusion>(baseModel);

    // Access the patch mobility, diffusivity and electric field
    const scalarField& muf = ddModel.mobility().muPatch(p.index());
    const scalarField& Df  = ddModel.diffusivity().DPatch(p.index());
    const fvPatchField<vector>& Ef = 
                                p.lookupPatchField<volVectorField, vector>("E");

    // Physics Calculations
    word scheme = ddModel.fluxScheme();
    const scalarField uDrift_n = Z * (muf * Ef) & nf;
    const tmp<scalarField> tUEff = this->calcEffectiveWallVelocity
                                                               (m, T, uDrift_n);
    const tmp<scalarField> tUAbs = this->calcAbsorptionVelocity(m, T, uDrift_n);
    const scalarField& uEff = tUEff();
    const scalarField& uAbs = tUAbs();

    // Set mixed b.c. parameters
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->operator==(this->patchInternalField());
    scalarField& f = this->valueFraction();

    // Set the D/Δ 
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
            const scalar Pe = uDrift_n[faceI] / (D_delta[faceI] + VSMALL);
            const scalar num = uEff[faceI];
            const scalar den = D_delta[faceI] * Bern(Pe) + uAbs[faceI];
            f[faceI] = num / (den + VSMALL);
        }
    }
    else // Standard schemes (e.g. upwind(div) + central(laplacian))
    {
        forAll(p, faceI)
        {
            f[faceI] = uEff[faceI] / (D_delta[faceI] + uEff[faceI] + SMALL);
        }
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}

void ddWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    // Write standard Mixed BC entries (valueFraction, refValue, etc.)
    mixedFvPatchScalarField::write(os);

    // Write our custom entries so the simulation can be restarted
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
