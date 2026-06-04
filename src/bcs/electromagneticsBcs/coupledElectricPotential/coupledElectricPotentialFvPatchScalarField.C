/*---------------------------------------------------------------------------*\
  File: coupledElectricPotentialFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::coupledElectricPotentialFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "IOField.H"
#include "mappedPatchFieldBase.H"

#include "coupledElectricPotentialFvPatchScalarField.H"
#include "SoPLASMAConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField&
coupledElectricPotentialFvPatchScalarField::getOrCreateField
(
    const word& fieldName
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    auto* ptr = mesh.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}

void coupledElectricPotentialFvPatchScalarField::storeECFields
(
    const word& prefix,
    const scalarField& sec,
    const scalarField& secPatch
)
const
{
    volScalarField& ec =
        getOrCreateField(IOobject::scopedName(prefix, "ec"));
    ec.boundaryFieldRef()[patch().index()] = sec;

    volScalarField& ecPatch =
        getOrCreateField(IOobject::scopedName(prefix, "ecPatch"));
    ecPatch.boundaryFieldRef()[patch().index()] = secPatch;
}

void coupledElectricPotentialFvPatchScalarField::writeFileHeader
(
    Ostream& os
)
{
    writeCommented(os, "Time");
    writeTabbed(os, "Q_[C]");
    writeTabbed(os, "J_[C/m^2]");
    writeTabbed(os, "Ec_[F/m^2]");
    writeTabbed(os, "EcPatch_[F/m^2]");
    writeTabbed(os, "phiMin_[V]");
    writeTabbed(os, "phiMax_[V]");
    writeTabbed(os, "phiAvg_[V]");
    writeTabbed(os, "phiNbrMin_[V]");
    writeTabbed(os, "phiNbrMax_[V]");
    writeTabbed(os, "phiNbrAvg_[V]");

    os << endl;

    writtenHeader_ = true;
    updateHeader_  = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),
    functionObjects::writeFile
    (
        db(),
        "coupledElectricPotential",
        "undefined",
        false
    ),
    phiNbrName_("undefined-phiNbr"),
    surfChargeNbrName_("undefined-surfChargeNbr"),
    surfChargeName_("undefined-surfCharge"),
    logInterval_(-1),
    executionIndex_(0),
    verbose_(false),
    prefix_()
{
    this->refValue()      = Zero;
    this->refGrad()       = Zero;
    this->valueFraction() = 1.0;
    this->source()        = 0.0;
}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    phiNbrName_(psf.phiNbrName_),
    surfChargeNbrName_(psf.surfChargeNbrName_),
    surfChargeName_(psf.surfChargeName_),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    functionObjects::writeFile
    (
        db(),
        "coupledElectricPotential",
        patch().name(),
        false
    ),
    phiNbrName_(dict.getOrDefault<word>("phiNbr", "ePotential")),
    surfChargeNbrName_(dict.getOrDefault<word>("surfChargeNbr", "none")),
    surfChargeName_(dict.getOrDefault<word>("surfCharge", "none")),
    logInterval_(dict.getOrDefault<scalar>("logInterval", -1)),
    executionIndex_(0),
    verbose_(dict.getOrDefault<bool>("verbose", false)),
    prefix_(dict.getOrDefault<word>("prefix", "multiWorld"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        refValue()      = *this;
        refGrad()       = Zero;
        valueFraction() = 1.0;
    }

    bool boolVal(false);
    if (dict.readIfPresent("useImplicit", boolVal))
    {
        this->useImplicit(boolVal);
    }

    if (dict.found("source"))
    {
        source() = scalarField("source", dict, p.size());
    }
    else
    {
        source() = 0.0;
    }

    writeFile::read(dict);
}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    phiNbrName_(psf.phiNbrName_),
    surfChargeNbrName_(psf.surfChargeNbrName_),
    surfChargeName_(psf.surfChargeName_),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf
)
:
    mixedFvPatchScalarField(psf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), psf.internalField()),
        *this,
        psf
    ),
    functionObjects::writeFile(psf),
    phiNbrName_(psf.phiNbrName_),
    surfChargeNbrName_(psf.surfChargeNbrName_),
    surfChargeName_(psf.surfChargeName_),
    logInterval_(psf.logInterval_),
    executionIndex_(psf.executionIndex_),
    verbose_(psf.verbose_),
    prefix_(psf.prefix_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledElectricPotentialFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
}


void coupledElectricPotentialFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


scalar coupledElectricPotentialFvPatchScalarField::epsilon() const
{
    scalar epsilonR = 1.0;
    
    if (db().foundObject<IOdictionary>("dielectricProperties"))
    {
        epsilonR = db().lookupObject<IOdictionary>("dielectricProperties")
            .getOrDefault<scalar>("dielectricConstant", 1.0);
    }
    else if (db().time().foundObject<electromagneticsModel>
    (
        "electromagneticsProperties"
    )
    )
    {
        epsilonR = db().time().lookupObject<electromagneticsModel>
        (
            "electromagneticsProperties"
        ).epsilonR();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find permittivity for region '"
            << db().name() << "'." << nl
            << "Expected either an electromagneticsModel registered as "
            << "'electromagneticsProperties' (gas region) or an "
            << "'dielectricProperties' dictionary (dielectric region)." << nl
            << exit(FatalError);
    }

    return epsilonR * constant::plasma::epsilon0.value();
}

void coupledElectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const int oldTag = UPstream::incrMsgType();

    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    const scalarField phiC(patchInternalField());
    const scalarField& phiP = *this;

    const scalar epsilonPhiP = epsilon();
    const scalarField epsilonDelta(epsilonPhiP * patch().deltaCoeffs());

    scalarField phiCNbr;
    scalarField phiPNbr;
    scalarField epsilonDeltaNbr;

    if (mpp.sameWorld())
    {
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const label samplePatchi = mpp.samplePolyPatch().index();
        const fvPatch& nbrPatch =
            refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

        const auto& nbrField =
            refCast<const coupledElectricPotentialFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField>(phiNbrName_)
            );

        phiCNbr = nbrField.patchInternalField();
        phiPNbr = nbrField;
        const scalar epsilonNbr = nbrField.epsilon();
        epsilonDeltaNbr = epsilonNbr * nbrPatch.deltaCoeffs();
    }
    else
    {
        phiCNbr = patchInternalField();
        phiPNbr = phiP;
        epsilonDeltaNbr = epsilonDelta;
    }

    distribute(this->internalField().name() + "_value", phiCNbr);
    distribute(this->internalField().name() + "_patchValue", phiPNbr);
    distribute(this->internalField().name() + "_weights", epsilonDeltaNbr);

    scalarField epsilonDeltaC(this->size(), GREAT);

    scalarField surfCharge(phiP.size(), Zero);
    if (surfChargeName_ != "none")
        {
            if (db().foundObject<volScalarField>(surfChargeName_))
            {
                const volScalarField& sigma = 
                    db().lookupObject<volScalarField>(surfChargeName_);
                
                surfCharge = sigma.boundaryField()[patch().index()];
            }
            else
            {
                WarningInFunction 
                    << "Field '" << surfChargeName_<< "' not found in registry." 
                    << endl;
            }
        }

    scalarField surfChargeNbr(phiP.size(), Zero);
    if (surfChargeNbrName_ != "none")
    {
        if (mpp.sameWorld())
        {
            const polyMesh& nbrMesh = mpp.sampleMesh();
            const label samplePatchi = mpp.samplePolyPatch().index();

            if (nbrMesh.foundObject<volScalarField>(surfChargeNbrName_))
            {
                const volScalarField& sigmaNbr = 
                    nbrMesh.lookupObject<volScalarField>
                    (surfChargeNbrName_);
                
                surfChargeNbr = sigmaNbr.boundaryField()[samplePatchi];
            }
        }
        distribute(surfChargeNbrName_, surfChargeNbr);
    }

    if (surfChargeName_ != "none" && surfChargeNbrName_ != "none")
    {
        if (gMax(mag(surfCharge)) > VSMALL && gMax(mag(surfChargeNbr)) > VSMALL)
        {
            WarningInFunction
                << "Both surfCharge and surfChargeNbr are non-zero "
                << "on this interface. This will add them and may "
                << "double-count surface charge.\n";
        }
    }
    
    valueFraction() = epsilonDeltaNbr/(epsilonDeltaNbr + epsilonDelta);
    refValue() = phiCNbr;
    refGrad() = (surfCharge + surfChargeNbr)/epsilonPhiP;
    source()  = Zero;

    if (this->useImplicit())
    {
        source() =
            epsilonDelta*patch().magSf() *
            (
                valueFraction()*deltaPhi()
            + (surfCharge + surfChargeNbr)/beta()
            );
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (verbose_)
    {
        // Total transferred charge across interface [C]
        const scalar Q = gSum(epsilonPhiP * patch().magSf() * snGrad());

        // Total interface area [m^2]
        const scalar magSf = gSum(patch().magSf());

        // Area-averaged normal electric flux density [C/m^2]
        const scalar q = Q / max(magSf, VSMALL);

        // Effective dielectric coupling [F/m^2] from Δφ (owner vs cell-center)
        const scalarField qField(epsilonPhiP * snGrad());
        const scalarField deltaPhi(phiCNbr - phiC);
        scalarField ec(deltaPhi.size(), Zero);

        forAll(deltaPhi, i)
        {
            if (mag(deltaPhi[i]) > VSMALL)
            {
                ec[i] = qField[i] / deltaPhi[i];
            }
        }
        const scalar aveEc =
            gSum(ec * patch().magSf()) / max(magSf, VSMALL);

        // Same coupling estimate from Δφ (patch values)
        const scalarField deltaPhiPatch(phiPNbr - phiP);
        scalarField ecPatch(deltaPhiPatch.size(), Zero);

        forAll(deltaPhiPatch, i)
        {
            if (mag(deltaPhiPatch[i]) > VSMALL)
            {
                ecPatch[i] = qField[i] / deltaPhiPatch[i];
            }
        }
        const scalar aveEcPatch =
            gSum(ecPatch * patch().magSf()) / max(magSf, VSMALL);

        // Potential statistics [V]
        const scalarMinMax phiPMinMax = gMinMax(phiP);
        const scalar phiPAvg = gAverage(phiP);
        const scalarMinMax phiPNbrMinMax = gMinMax(phiPNbr);
        const scalar phiPNbrAvr = gAverage(phiPNbr);

        Info<< nl
            << patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :" << nl
            << " Total transferred charge [C]:" << Q << nl
            << " Area [m^2]:" << magSf << nl
            << " Normal electric flux density [C/m^2]:" << q << nl
            << " Effective dielectric coupling [F/m^2] (cell-center):"
            << aveEc << nl
            << " Effective dielectric coupling [F/m^2] (patch-based):"
            << aveEcPatch << nl
            << " Local potential [V]"
            << " min:" << phiPMinMax.min()
            << " max:" << phiPMinMax.max()
            << " avg:" << phiPAvg << nl
            << " Neighbour potential [V]"
            << " min:" << phiPNbrMinMax.min()
            << " max:" << phiPNbrMinMax.max()
            << " avg:" << phiPNbrAvr
            << nl << endl;

       // Handle data for file output
        if (canResetFile())
        {
            resetFile(patch().name());
        }

        if (canWriteHeader())
        {
            writeFileHeader(file());
        }

        if (canWriteToFile() && writeFile())
        {
            file()
                << db().time().timeOutputValue() << token::TAB
                << Q << token::TAB
                << q << token::TAB
                << aveEc << token::TAB
                << aveEcPatch << token::TAB
                << phiPMinMax.min() << token::TAB
                << phiPMinMax.max() << token::TAB
                << phiPAvg << token::TAB
                << phiPNbrMinMax.min() << token::TAB
                << phiPNbrMinMax.max() << token::TAB
                << phiPNbrAvr << token::TAB
                << endl;
        }


        // Store etc fields as patch fields of a volScalarField
        storeECFields(prefix_, ec, ecPatch);
    }

    UPstream::msgType(oldTag);  // Restore tag
}

void coupledElectricPotentialFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const label mat,
    const direction cmpt
)
{
    label index = this->patch().index();

    const label nbrPatchId =  this->patch().patch().neighbPolyPatchID();

    const label globalPatchID =
        matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

    const label meshNrbId = matrix.lduMeshAssembly().findNbrMeshId
    (
        this->patch().patch(),
        mat
    );

    const mixedFvPatchScalarField& fPatch =
        static_cast<const mixedFvPatchScalarField&>(*this);

    const Field<scalar> intCoeffsCmpt
    (
        matrix.internalCoeffs()[globalPatchID].component(cmpt)
    );

    const scalarField sourceCorr(fPatch.source());

    const labelList& faceMap =
        matrix.lduMeshAssembly().faceBoundMap()[mat][index];

    const labelList& myCells =
        matrix.lduMeshAssembly().cellBoundMap()[meshNrbId][nbrPatchId];

    const labelList& nbrCells =
        matrix.lduMeshAssembly().cellBoundMap()[mat][index];

    const word patchName = this->patch().name();

    forAll(faceMap, j)
    {
        label globalFaceI = faceMap[j];
        label myCellI     = myCells[j];
        label nbrCellI    = nbrCells[j];

        const scalar intCorr = -intCoeffsCmpt[j];
        const scalar srcCorr = -sourceCorr[j];

        if (this->patch().patch().masterImplicit())
        {
            if (myCellI > nbrCellI)
            {
                if (matrix.asymmetric())
                {
                    matrix.lower()[globalFaceI] += intCorr;
                }
            }
            else
            {
                matrix.upper()[globalFaceI] += intCorr;

            }

            matrix.diag()[myCellI] -= intCorr;

            matrix.source()[myCellI] += srcCorr;
        }
        else
        {
            if (myCellI < nbrCellI)
            {
                matrix.upper()[globalFaceI] += intCorr;
            }
            else
            {
                if (matrix.asymmetric())
                {
                    matrix.lower()[globalFaceI] += intCorr;
                }
            }

            matrix.diag()[myCellI] -= intCorr;

            matrix.source()[myCellI] += srcCorr;
        }
    }
}

void coupledElectricPotentialFvPatchScalarField::initEvaluate
(
    const UPstream::commsTypes commsType
)
{

    this->setUpdated(false);

    mixedFvPatchScalarField::initEvaluate(commsType);
}

void coupledElectricPotentialFvPatchScalarField::evaluate
(
    const UPstream::commsTypes commsType
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    scalarField own  = patchInternalField();
    scalarField temp = own + refGrad()/patch().deltaCoeffs();

    scalarField newValues =
        valueFraction()*refValue()
        + (scalar(1) - valueFraction()) * temp;

    Field<scalar>::operator=(newValues);

    fvPatchField<scalar>::evaluate(commsType);
}

tmp<scalarField> coupledElectricPotentialFvPatchScalarField::
beta() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    if (!mpp.sameWorld())
    {
        FatalErrorInFunction
            << "Coupled electric potential not supported in combination "
            << "with multi-world setups."
            << exit(FatalError);
    }

    const label samplePatchi = mpp.samplePolyPatch().index();
    const polyMesh& nbrMesh = mpp.sampleMesh();

    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    const coupledElectricPotentialFvPatchScalarField&
        nbrField = refCast
            <const coupledElectricPotentialFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField>(phiNbrName_)
            );

    scalarField phiCNbr(nbrField.patchInternalField());
    mpp.distribute(phiCNbr);

    scalarField epsilonDeltaNbr(nbrField.epsilon() * nbrPatch.deltaCoeffs());
    mpp.distribute(epsilonDeltaNbr);

    scalarField epsilonDelta(epsilon() * patch().deltaCoeffs());

    return (epsilonDeltaNbr + epsilonDelta);
}

tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::deltaPhi() const
{
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    if (!mpp.sameWorld())
    {
        FatalErrorInFunction
            << "Coupled electric potential not supported in combination "
            << "with multi-world setups."
            << exit(FatalError);
    }

    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    const coupledElectricPotentialFvPatchScalarField& nbrField =
        refCast<const coupledElectricPotentialFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField>(phiNbrName_)
        );

    const scalarField& phiLocal = *this;
    scalarField phiNbr(nbrField);
    mpp.distribute(phiNbr);

    return tmp<scalarField>::New(phiNbr - phiLocal);
}


bool coupledElectricPotentialFvPatchScalarField::writeFile()
{
    if (!verbose_ || (logInterval_ <= 0))
    {
        return false;
    }

    const auto& time = patch().boundaryMesh().mesh().time();

    const scalar t = time.timeOutputValue();
    const scalar ts = time.startTime().value();
    const scalar deltaT = time.deltaTValue();

    const label executionIndex = label
    (
        (
            (t - ts)
          + 0.5*deltaT
        )
        /logInterval_
    );

    bool write = false;
    if (executionIndex > executionIndex_)
    {
        executionIndex_ = executionIndex;
        write = true;
    }

    return write;
}


void coupledElectricPotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);

    os.writeEntryIfDifferent<word>("phiNbr", "ePotential", phiNbrName_);
    os.writeEntryIfDifferent<word>("surfChargeNbr", "none", surfChargeNbrName_);
    os.writeEntryIfDifferent<word>("surfCharge", "none", surfChargeName_);
    os.writeEntry<scalar>("logInterval", logInterval_);

    os.writeEntryIfDifferent<bool>("verbose", false, verbose_);
    os.writeEntryIfDifferent<word>("prefix", "multiWorld", prefix_);

    // Write writeFile entries
    os.writeEntry<label>("writePrecision", writePrecision_);
    os.writeEntry<bool>("updateHeader", updateHeader_);
    os.writeEntry<bool>("writeToFile", writeToFile_);
    os.writeEntry<bool>("useUserTime", useUserTime_);

    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    coupledElectricPotentialFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
