/*---------------------------------------------------------------------------*\
  File: electronDDWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electronDDWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "plasmaTransport.H"
#include "electronDDWallFluxImplicitFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    electronDDWallFluxImplicitFvPatchScalarField
);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void electronDDWallFluxImplicitFvPatchScalarField::buildSEECList() const
{
    if (mapped_) return;

    if (!db().foundObject<plasmaTransport>("plasmaTransport"))
    {
        return;
    }

    const plasmaTransport& transport =
                          db().lookupObject<plasmaTransport>("plasmaTransport");

    const plasmaSpecies& speciesDB = transport.species();

    seec_.setSize(speciesDB.nSpecies(), defaultSEEC_);

    forAll(speciesDB.speciesNames(), specI)
    {
        const word& name = speciesDB.speciesNames()[specI];
        if (speciesSEEC_.found(name))
        {
            seec_[specI] =
                speciesSEEC_.get<scalar>(name);
        }
    }

    mapped_ = true;
}

tmp<scalarField>
electronDDWallFluxImplicitFvPatchScalarField::calcSEEFlux() const
{
    const fvPatch& p = patch();

    // Ensure the SEEC list is built
    buildSEECList();

    if (!mapped_)
    {
        return tmp<scalarField>::New(p.size(), Zero);
    }

    const plasmaTransport& transport =
                          db().lookupObject<plasmaTransport>("plasmaTransport");
    const plasmaSpecies& speciesDB = transport.species();
    const scalarField& magSf = p.magSf();

    tmp<scalarField> tSEE = tmp<scalarField>::New(p.size(), Zero);
    scalarField& SEE = tSEE.ref();

    for (const label specI : speciesDB.positiveIonSpeciesIDs())
    {
        const word fluxName =
            "particleFlux_" + speciesDB.speciesName(specI);

        if (!p.boundaryMesh().mesh().foundObject<surfaceScalarField>(fluxName))
            continue;

        const fvsPatchScalarField& phiI =
            p.lookupPatchField<surfaceScalarField, scalar>(fluxName);

        SEE += seec_[specI] * max(0.0, phiI / magSf);
    }

    return tSEE;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
electronDDWallFluxImplicitFvPatchScalarField::calcWallThermalVelocity
(
    const dimensionedScalar& m,
    const scalarField& T
) const
{
    // Pure thermal velocity u_th/4 only.
    // Drift correction is handled separately in valueInternalCoeffs().
    return calcThermalVelocity(m, T);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
electronDDWallFluxImplicitFvPatchScalarField::
electronDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF),
    enableSurfaceCharging_(true),
    includeDriftFlux_(false),
    enableSEE_(false),
    defaultSEEC_(0.0),
    speciesSEEC_(dictionary::null),
    seec_(0),
    mapped_(false)
{}

// Dictionary Constructor
electronDDWallFluxImplicitFvPatchScalarField::
electronDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_(dict.lookupOrDefault<bool>
        ("enableSurfaceCharging", true)),
    includeDriftFlux_(dict.lookupOrDefault<bool>("includeDriftFlux", false)),
    enableSEE_(dict.lookupOrDefault<bool>("enableSEE", false)),
    defaultSEEC_(dict.lookupOrDefault<scalar>("defaultSEEC", 0.05)),
    speciesSEEC_(dict.subOrEmptyDict("speciesSEEC")),
    seec_(0),
    mapped_(false)
{}

// Mapping Constructor
electronDDWallFluxImplicitFvPatchScalarField::
electronDDWallFluxImplicitFvPatchScalarField
(
    const electronDDWallFluxImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor
electronDDWallFluxImplicitFvPatchScalarField::
electronDDWallFluxImplicitFvPatchScalarField
(
    const electronDDWallFluxImplicitFvPatchScalarField& ptf
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor (with new internal field)
electronDDWallFluxImplicitFvPatchScalarField::
electronDDWallFluxImplicitFvPatchScalarField
(
    const electronDDWallFluxImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<scalar>>
electronDDWallFluxImplicitFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& weights
) const
{
    // If drift is disabled: return 0 always
    // fvm::div contributes nothing. All flux via fvm::laplacian thermal term
    if (!includeDriftFlux_)
    {
        return tmp<Field<scalar>>
        (
            new scalarField(this->size(), Zero)
        );
    }

    // If drift is enabled: delegate to base class
    // Base class returns neg(Z*mu*E_n): 1 if toward wall, 0 if away
    return ddWallFluxImplicitFvPatchScalarField::valueInternalCoeffs(weights);
}

tmp<Field<scalar>>
electronDDWallFluxImplicitFvPatchScalarField::gradientInternalCoeffs() const
{
    // Thermal flux -(u_th/4)/D. Always present, same as base class
    return ddWallFluxImplicitFvPatchScalarField::gradientInternalCoeffs();
}

tmp<Field<scalar>>
electronDDWallFluxImplicitFvPatchScalarField::gradientBoundaryCoeffs() const
{
    // No SEE: return zero (no explicit RHS source)
    if (!enableSEE_)
    {
        return tmp<Field<scalar>>
        (
            new scalarField(this->size(), Zero)
        );
    }

    const fvPatch& p = patch();

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

    // SEE flux [particles/m^2/s]
    tmp<scalarField> tSEE = calcSEEFlux();
    const scalarField& SEEflux = tSEE();

    // SEE source: gradBoundaryCoeff = -SEEflux / D
    // (negative so laplacian adds positive source contribution)
    return tmp<Field<scalar>>
    (
        new scalarField(-SEEflux / (Df + VSMALL))
    );
}

void electronDDWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxImplicitFvPatchScalarField::write(os);   

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
    os.writeEntry("includeDriftFlux", includeDriftFlux_);
    os.writeEntry("enableSEE", enableSEE_);
    os.writeEntry("defaultSEEC", defaultSEEC_);
    if (!speciesSEEC_.empty())
    {
        os.writeEntry("speciesSEEC", speciesSEEC_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
