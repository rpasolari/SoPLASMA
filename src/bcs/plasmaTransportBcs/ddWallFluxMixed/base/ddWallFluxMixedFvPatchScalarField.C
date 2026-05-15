/*---------------------------------------------------------------------------*\
  File: ddWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ddWallFluxMixedFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvPatchFieldMapper.H"

#include "plasmaTransport.H"
#include "plasmaSpecies.H"
#include "ddWallFluxMixedFvPatchScalarField.H"
#include "SoPLASMAConstants.H"
#include "driftDiffusion.H"

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

// void ddWallFluxMixedFvPatchScalarField::updateCoeffs()
// {
//     if (updated()) return;

//     const fvPatch& p = patch();
//     if (p.size() == 0)
//     {
//         mixedFvPatchScalarField::updateCoeffs();
//         return;
//     }
//     const vectorField nf = p.nf();
//     const scalarField& delta = p.deltaCoeffs();

//     // Determine the species name from the field itself (e.g., n_e -> e)
//     const word fieldName = this->internalField().name();
//     word speciesName = fieldName;
//     if (speciesName.startsWith("n_"))
//     {
//         speciesName.erase(0, 2);
//     }

//     tmp<scalarField> tT;

//     if (TName_ == "none")
//     {
//         tT = tmp<scalarField>::New(p.size(), TValue_.value());
//     }
//     else if (db().foundObject<volScalarField>(TName_))
//     {
//         tT = p.lookupPatchField<volScalarField, scalar>(TName_);
//     }
//     else if (db().foundObject<dimensionedScalar>(TName_))
//     {
//         const auto& udf = db().lookupObject<UniformDimensionedField<scalar>>
//                                                                        (TName_);
//         tT = tmp<scalarField>::New(p.size(), udf.value());
//     }
//     else
//     {
//         FatalErrorInFunction
//             << "Temperature identifier '" << TName_ << "' was not found in the "
//             << "object registry as a volScalarField or a dimensionedScalar."
//             << exit(FatalError);
//     }

//     const scalarField& T = tT();

//     // Patch based lookups for mu, D, and E
//     const fvPatchField<scalar>& muf = 
//         p.lookupPatchField<volScalarField, scalar>("mu_" + speciesName);

//     const fvPatchField<scalar>& Df = 
//         p.lookupPatchField<volScalarField, scalar>("D_" + speciesName);

//     const fvPatchField<vector>& Ef = 
//         p.lookupPatchField<volVectorField, vector>("E");

//     // Registry Lookups for Mass and Scheme
//     const plasmaTransport& transport = 
//         db().lookupObject<plasmaTransport>("plasmaTransport");
//     const plasmaSpecies& speciesDB = transport.species();
//     const label speciesID = speciesDB.speciesID(speciesName);
//     const dimensionedScalar& m = speciesDB.speciesMass(speciesID);
//     const scalar& Z = speciesDB.speciesChargeNumber(speciesID);

//     const plasmaTransportModel& model = transport.model(speciesID);
//     word scheme = "standard";
//     if (isA<driftDiffusion>(model))
//     {
//         const driftDiffusion& ddModel = refCast<const driftDiffusion>(model);
//         scheme = ddModel.fluxScheme();
//     }

//     // Calculate Normal Drift Velocity (mu * E_normal)
//     // Note: Child classes handle the charge sign inside calcWallVelocity
//     scalarField uDrift_n = (muf * Ef) & nf;

//     // Get the wall velocity from child
//     tmp<scalarField> tUWall = this->calcWallVelocity(m, T, uDrift_n, Z);
//     const scalarField& uWall = tUWall();

//     // Reset standard Mixed BC parameters
//     this->refValue() = 0.0;
//     this->refGrad() = 0.0;
//     this->operator==(this->patchInternalField());
//     scalarField& f = this->valueFraction();

//     // Define the D / distance parameter
//     scalarField D_delta = Df * delta;

//     if (scheme == "ScharfetterGummel")
//     {
//         // Bernoulli Function Definition: B(x) = x / exp(x) - 1
//         auto Bern = [](scalar x) -> scalar
//         {
//             const scalar ax = mag(x);
//             if (ax < 1e-4)
//             {
//                 return 1.0 - 0.5*x + (x*x)/12.0 - pow4(x)/720.0;;
//             }
//             if (x > 200.0)  return 0.0;
//             if (x < -200.0) return -x;
//             return x / (Foam::exp(x) - 1.0);
//         };

//         forAll(p, faceI)
//         {
//             scalar Pe = uDrift_n[faceI] / (D_delta[faceI] + VSMALL);

//             // Use of the identity Be(Pe) - Be(-Pe) = -Pe
//             scalar num = uWall[faceI] - D_delta[faceI] * Pe;
//             scalar den = D_delta[faceI] * Bern(Pe) + uWall[faceI];

//             f[faceI] = num / (den + VSMALL);
//         }
//     }
//     else    //standard schemes
//     {
//         f = uWall / (D_delta + uWall + VSMALL);
//     }

//     mixedFvPatchField<scalar>::updateCoeffs();
// }

void ddWallFluxMixedFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;

    const fvPatch& p = patch();
    if (p.size() == 0)
    {
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    const fvMesh& mesh = p.boundaryMesh().mesh(); // Access the mesh from the patch
    const vectorField nf = p.nf();
    const scalarField& delta = p.deltaCoeffs();

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

void ddWallFluxMixedFvPatchScalarField::write(Ostream& os) const
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
