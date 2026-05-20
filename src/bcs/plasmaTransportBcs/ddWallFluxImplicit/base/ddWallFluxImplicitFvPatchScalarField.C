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

#include "ddWallFluxImplicitFvPatchScalarField.H"

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
    fvPatchScalarField(p, iF),
    TName_("none"),
    TValue_("T", dimTemperature, 300.0)
{}

// Dictionary Constructor
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchScalarField(p, iF, dict),
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
    fvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// Copy Constructor (from another patch field)
ddWallFluxImplicitFvPatchScalarField::ddWallFluxImplicitFvPatchScalarField
(
    const ddWallFluxImplicitFvPatchScalarField& ptf
)
:
    fvPatchScalarField(ptf),
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
    fvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    TValue_(ptf.TValue_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ddWallFluxImplicitFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Set n_face = n_cell (consistent with valueInternalCoeffs = 1 for drift)
    // The actual flux is controlled by the coefficient methods, not n_face.
    // This ensures post-processing and fvc operators see a sensible face value.
    fvPatchScalarField::operator==(this->patchInternalField());

    fvPatchField<scalar>::evaluate();
}

void ddWallFluxImplicitFvPatchScalarField::updateCoeffs()
{
    // Coefficients are computed on-the-fly in the coefficient methods.
    // Nothing to precompute here.
    if (updated()) return;
    fvPatchScalarField::updateCoeffs();
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvPatch& p = patch();

    // Access the normal vector and delta coeffs
    const vectorField nf = p.nf();

    // Determine species name (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

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

    // Access the patch mobility and electric field
    const scalarField& muf = ddModel.mobility().muPatch(p.index());
    const fvPatchField<vector>& Ef = 
                                p.lookupPatchField<volVectorField, vector>("E");

    const scalarField uDrift_n = Z * (muf * Ef) & nf;

    // Implicit drift flux: upwind coeff = 1 if drift toward wall, else 0
    return tmp<Field<scalar>>
    (
        new scalarField(pos(uDrift_n))
    );
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // No explicit RHS contribution from drift
    return tmp<Field<scalar>>
    (
        new scalarField(this->size(), Zero)
    );
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::gradientInternalCoeffs() const
{
    const fvPatch& p = patch();

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

    // Access the patch diffusivity
    const scalarField& Df  = ddModel.diffusivity().DPatch(p.index());

    // Physics Calculations
    tmp<scalarField> tUth = this->calcWallThermalVelocity(m, T);
    const scalarField& uth4 = tUth();

    // Implicit thermal flux: gradCoeff = -u_th/4 / D
    // (negative so fvm::laplacian adds positive diagonal contribution)
    return tmp<Field<scalar>>
    (
        new scalarField(-uth4 / (Df + VSMALL))
    );
}

tmp<Field<scalar>>
ddWallFluxImplicitFvPatchScalarField::gradientBoundaryCoeffs() const
{
    // No explicit RHS contribution from diffusion
    return tmp<Field<scalar>>
    (
        new scalarField(this->size(), Zero)
    );
}

void ddWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    // Write our custom entries so the simulation can be restarted
    if (TName_ != "none")
    {
        os.writeEntry("T", TName_);
    }
    else
    {
        os.writeEntry("T", TValue_.value());
    }

    fvPatchScalarField::writeValueEntry(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
