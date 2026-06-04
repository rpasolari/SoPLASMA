/*---------------------------------------------------------------------------*\
  File: electronDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electronDDWallFluxMixedFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "plasmaTransport.H"
#include "electronDDWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    electronDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> 
electronDDWallFluxMixedFvPatchScalarField::calcAbsorptionVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    tmp<scalarField> tVel = calcThermalVelocity(m, T);
    scalarField& uAbs = tVel.ref();

    // If drift flux is enabled, add the directed motion component
    if (includeDriftFlux_)
    {
        uAbs += max(0.0, uDriftNormal);
    }
    
    return tVel;
}

tmp<scalarField> 
electronDDWallFluxMixedFvPatchScalarField::calcEffectiveWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal
) const
{
    tmp<scalarField> tVel = calcThermalVelocity(m, T);
    scalarField& uWall = tVel.ref();

    // If drift flux is enabled, add the directed motion component
    if (includeDriftFlux_)
    {
        uWall += max(0.0, -uDriftNormal);
    }
    else
    {
        uWall -= uDriftNormal;
    }

    return tVel;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF),
    enableSurfaceCharging_(false),
    includeDriftFlux_(false),
    enableSEE_(false),
    defaultSEEC_(0.0),
    speciesSEEC_(dictionary::null),
    seec_(0),
    mapped_(false)
{}

// Dictionary Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_
            (dict.lookupOrDefault<bool>("enableSurfaceCharging", false)),
    includeDriftFlux_(dict.lookupOrDefault<bool>("includeDriftFlux", false)),
    enableSEE_(dict.lookupOrDefault<bool>("enableSEE", false)),
    defaultSEEC_(dict.lookupOrDefault<scalar>("defaultSEEC", 0.05)),
    speciesSEEC_(dict.subOrEmptyDict("speciesSEEC")),
    seec_(0), 
    mapped_(false)
{}

// Mapping Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf
)
:
    ddWallFluxMixedFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// Copy Constructor (with new internal field)
electronDDWallFluxMixedFvPatchScalarField::
electronDDWallFluxMixedFvPatchScalarField
(
    const electronDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_),
    includeDriftFlux_(ptf.includeDriftFlux_),
    enableSEE_(ptf.enableSEE_),
    defaultSEEC_(ptf.defaultSEEC_),
    speciesSEEC_(ptf.speciesSEEC_),
    seec_(ptf.seec_),
    mapped_(ptf.mapped_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void electronDDWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Initialize the SEE entries
    if (enableSEE_ && !mapped_)
    {
        if (!db().foundObject<plasmaTransport>("plasmaTransport"))
        {
            ddWallFluxMixedFvPatchScalarField::updateCoeffs();
            return;
        }

        const plasmaSpecies& speciesDB =
            db().lookupObject<plasmaTransport>("plasmaTransport").species();

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

    // Call the standard updateCoeffs from the base class
    ddWallFluxMixedFvPatchScalarField::updateCoeffs();

    if (!enableSEE_) return;

    // Determine the species name from the field itself (e.g., n_e -> e)
    const word fieldName = this->internalField().name();
    word speciesName = fieldName;
    if (speciesName.startsWith("n_"))
    {
        speciesName.erase(0, 2);
    }

    // Modify refValue or refGrad for secondary electron emission
    const fvPatch& p = patch();

    const scalarField& magSf = p.magSf();
    scalarField totalSEE(p.size(), 0.0);

    // Registry lookup for the plasmaTransport object
    const plasmaTransport& transport = 
                        db().lookupObject<plasmaTransport>("plasmaTransport");

    const plasmaSpecies& speciesDB = transport.species();

    for (const label specI : speciesDB.positiveIonSpeciesIDs())
    {
        const word fluxName =
            "particleFlux_" + speciesDB.speciesName(specI);

        if (!p.boundaryMesh().mesh().foundObject<surfaceScalarField>(fluxName))
            continue;

        const fvsPatchScalarField& phiI =
            p.lookupPatchField<surfaceScalarField, scalar>(fluxName);

        totalSEE += seec_[specI] * max(0.0, phiI / magSf);
    }

    const fvPatchField<scalar>& Df = 
                 p.lookupPatchField<volScalarField, scalar>("D_" + speciesName);

    this->refGrad() += totalSEE / (Df + VSMALL);
}

void electronDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxMixedFvPatchScalarField::write(os);   

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
