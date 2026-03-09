/*---------------------------------------------------------------------------*\
  File: solidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::solidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "plasmaTransport.H"
#include "solidSurfaceFluxFvPatchScalarField.H"
#include "foamPlasmaToolkitConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(solidSurfaceFluxFvPatchScalarField, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

dimensionedScalar solidSurfaceFluxFvPatchScalarField::calcThermalVelocity
(
    const dimensionedScalar& m,
    const dimensionedScalar& T
) const
{
    return sqrt 
    (
        (8.0 * constant::plasma::kappaBoltzmann * T)
        /
        (constant::mathematical::pi * m)
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
solidSurfaceFluxFvPatchScalarField::solidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    T_("T", dimTemperature, 300.0),
    enableSurfaceCharging_(false)
{
    this->refValue()      = Zero;
    this->refGrad()       = Zero;
    this->valueFraction() = 0.0;
}

// Dictionary Constructor
solidSurfaceFluxFvPatchScalarField::solidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    T_
    (
        dict.lookupOrDefault
        (
            "T", 
            dimensionedScalar
            (
                "T",
                dimTemperature,
                300.0
            )
        )
    ),
    enableSurfaceCharging_
    (
        dict.lookupOrDefault<bool>
        (
            "enableSurfaceCharging", 
            false
        )
    )
{
    if (this->readMixedEntries(dict))
    {
        // Full restart or values provided in dictionary
    }
    else
    {
        this->refValue()      = Zero;
        this->refGrad()       = Zero;
        this->valueFraction() = 0.0;
    }
}

// Mapping Constructor
solidSurfaceFluxFvPatchScalarField::solidSurfaceFluxFvPatchScalarField
(
    const solidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    T_(ptf.T_),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor (from another patch field)
solidSurfaceFluxFvPatchScalarField::solidSurfaceFluxFvPatchScalarField
(
    const solidSurfaceFluxFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    T_(ptf.T_),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor (from patch field and new internal field)
solidSurfaceFluxFvPatchScalarField::solidSurfaceFluxFvPatchScalarField
(
    const solidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    T_(ptf.T_),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidSurfaceFluxFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const fvPatch& p = patch();
    if (p.size() == 0)
    {
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }
    const vectorField nf = p.nf();
    const scalarField& delta = p.deltaCoeffs();

    // Find plasma transport from registry
    const plasmaTransport& transport = 
        db().lookupObject<plasmaTransport>("plasmaTransport");

    // Find plasma species
    const plasmaSpecies& speciesDB = transport.species();
    const word fieldName = this->internalField().name();
    const word speciesName = fieldName.substr(2);
    const label speciesID = speciesDB.speciesID(speciesName);

    // Access model physics
    const plasmaTransportModel& model = transport.model(speciesID);
    const volVectorField& driftVelocity = model.driftVelocity();
    const volScalarField* Dptr = model.diffusivity();
    const word& scheme = model.fluxScheme();

    const dimensionedScalar& q = speciesDB.speciesCharge(speciesID);
    const dimensionedScalar& m = speciesDB.speciesMass(speciesID);

    // Get field values on the patch
    const fvPatchField<vector>& uDrift = 
        driftVelocity.boundaryField()[p.index()];

    // Normal drift velocity
    scalarField uDrift_n = uDrift & nf; 

    scalarField uWall = this->calcWallVelocity(m, T_, uDrift_n);

    // Reset standard Mixed BC parameters
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    // this->valueFraction() = 0.0;
    this->operator==(this->patchInternalField());
    scalarField& f = this->valueFraction();

    if (Dptr)
    {
        // D exists -> Drift-Diffusion model is used then
        const fvPatchField<scalar>& D = Dptr->boundaryField()[p.index()];

        // Diffusivity / Delta
        scalarField D_delta = D * delta;

        if (scheme == "ScharfetterGummel")
        {
            // Bernoulli Function Definition: B(x) = x / exp(x) - 1
            auto Be = [](scalar x) -> scalar
            {
                if (mag(x) < 1e-4) 
                {
                    return 1.0 - x/2.0 + sqr(x)/12.0 - pow4(x)/720.0;
                }
                if (x > 200)    return 0.0;
                if (x < -200)   return -x;
                
                return x / (exp(x) - 1.0); 
            };

            forAll(p, faceI)
            {
                scalar Pe = uDrift_n[faceI] / (D_delta[faceI] + VSMALL);

                scalar num = D_delta[faceI] * Be(-Pe);
                scalar den = D_delta[faceI] * Be(Pe) + uWall[faceI];

                scalar K = num / (den + VSMALL);

                f[faceI] = 1.0 - K;
            }
        }
        else
        {
            // D is null -> Full momentum model is used then
            scalarField flowIn = pos0(uDrift_n);
            scalarField flowOut = 1.0 - flowIn;

            scalarField denom = uWall + D_delta - (uDrift_n * flowOut);

            scalarField num = uWall - uDrift_n;

            this->valueFraction() = num / (denom + VSMALL);
        }
    }
    else
    {
        this->valueFraction() = 0.0;
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}

bool solidSurfaceFluxFvPatchScalarField::enableSurfaceCharging() const
{
    return enableSurfaceCharging_;
}

void solidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    // 1. Write standard Mixed BC entries (valueFraction, refValue, etc.)
    mixedFvPatchScalarField::write(os);

    // 2. Write our custom entries so the simulation can be restarted
    os.writeEntry("T", T_);
    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
