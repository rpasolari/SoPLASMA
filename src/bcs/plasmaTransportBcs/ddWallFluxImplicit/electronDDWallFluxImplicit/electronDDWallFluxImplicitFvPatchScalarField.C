/*---------------------------------------------------------------------------*\
  File: electronDDWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electronDDWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> electronDDWallFluxImplicitFvPatchScalarField::calcWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal,
    const scalar Z
) const
{
    tmp<scalarField> tVel = calcThermalVelocity(m, T);
    scalarField& uWall = tVel.ref();

    // If drift flux is enabled, add the directed motion component
    if (includeDriftFlux_)
    {
        uWall += max(0.0, uDriftNormal);
    }
    else
    {
        uWall += uDriftNormal;
    }

    return tVel;
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
    setRefValue_(false),
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
    setRefValue_(dict.lookupOrDefault<bool>("setRefValue", false)),
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
    setRefValue_(ptf.setRefValue_),
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
    setRefValue_(ptf.setRefValue_),
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
    setRefValue_(ptf.setRefValue_),
    mapped_(ptf.mapped_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void electronDDWallFluxImplicitFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Initialize the SEE entries
    if (enableSEE_ && !mapped_)
    {
        if (db().foundObject<plasmaSpecies>("plasmaSpecies"))
        {
            const plasmaSpecies& speciesDB = 
                db().lookupObject<plasmaSpecies>("plasmaSpecies");

            seec_.setSize(speciesDB.nSpecies(), defaultSEEC_);

            forAll(speciesDB.speciesNames(), specI)
            {
                const word& name = speciesDB.speciesNames()[specI];
                
                if (speciesSEEC_.found(name))
                {
                    seec_[specI] = speciesSEEC_.lookupOrDefault<scalar>
                                        (name, defaultSEEC_);
                }
            }

            mapped_ = true;
        }
        else
        {
            // Run this if the plasmaSpecies registry is not ready
            ddWallFluxImplicitFvPatchScalarField::updateCoeffs();
            return;
        }
    }

    // Call the standard updateCoeffs from the base class
    ddWallFluxImplicitFvPatchScalarField::updateCoeffs();

    // Determine the species name from the field itself (e.g., n_e -> e)
    if (enableSEE_)
    {
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

        const wordList& names = speciesDB.speciesNames();

        for (label specI = 0; specI < names.size(); ++specI)
        {
            if (speciesDB.speciesChargeNumber(specI) > 0)
            {
                const word& ionName = speciesDB.speciesNames()[specI];
                const word fluxName = "particleFlux_" + ionName;

                const fvsPatchScalarField& phiI = 
                    p.lookupPatchField<surfaceScalarField, scalar>(fluxName);

                totalSEE += seec_[specI] * max(0.0, phiI / magSf);
            }
        }

        if (!setRefValue_)
        {
            const fvPatchField<scalar>& Df = 
                p.lookupPatchField<volScalarField, scalar>("D_" + speciesName);

            this->refGrad() += totalSEE / (Df + VSMALL);
        }

        this->evaluate();
    }
}

void electronDDWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxImplicitFvPatchScalarField::write(os);   

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
    os.writeEntry("includeDriftFlux", includeDriftFlux_);
    os.writeEntry("enableSEE", enableSEE_);
    os.writeEntry("defaultSEEC_", defaultSEEC_);
    if (!speciesSEEC_.empty())
    {
        os.writeEntry("speciesSEEC_", speciesSEEC_);
    }
    os.writeEntry("setRefValue_", setRefValue_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
