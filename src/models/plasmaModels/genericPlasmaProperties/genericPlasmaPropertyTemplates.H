/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ConstantProperty
    Foam::TabulatedProperty1D
    Foam::TabulatedProperty2D
    Foam::PowerLawProperty
    Foam::Function1Property

Description:
    Generic templated implementations for plasma physical properties.

    This header provides a library of mathematical logic that is independent of
    the specific physical quantity being modeled (e.g., Mobility, Diffusivity,
    Thermal Conductivity etc.).

    By using template specialization, these models are "stamped" onto physical 
    base classes (like plasmaMobilityModel). This ensures code reuse, 
    mathematical consistency across different physics, and centralized 
    maintenance.

\*---------------------------------------------------------------------------*/

#ifndef genericPlasmaPropertyTemplates_H
#define genericPlasmaPropertyTemplates_H

#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "Function1.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ConstantProperty Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class ConstantProperty
:
    public BaseModel
{
    scalar val_;

public:
    ConstantProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    )
    :
        BaseModel(modelType, dict, mesh, species, specieIndex),
        val_(0.0)
    {
        const word coeffKey = BaseModel::coeffKeyword();
        const word& sName = this->species_.speciesNames()[this->specieIndex_];

        if (!dict.found(coeffKey))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " 
                << BaseModel::modelType() << " model '" << modelType
                << "' requires entry '" << coeffKey << "'." << nl
                << exit(FatalIOError);
        }

        dict.readEntry(coeffKey, val_);

        this->checkValue(coeffKey, val_);
    }

    virtual void correct(volScalarField& field) const override
    {
        field = dimensionedScalar
                    (BaseModel::coeffKeyword(), field.dimensions(), val_);
        field.correctBoundaryConditions();
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty1D Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class TabulatedProperty1D
:
    public BaseModel
{
    //- The interpolation table (handles file reading and interpolation)
    interpolationTable<scalar> table_;

    //- The name of the field used as the lookup coordinate (e.g., "E" or "T_e")
    word lookupVariable_;

public:

    //- Constructor
    TabulatedProperty1D
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    )
    :
        BaseModel(modelType, dict, mesh, species, specieIndex),
        table_(dict),
        lookupVariable_(dict.get<word>("lookupVariable"))
    {
        const word& sName = this->species_.speciesNames()[this->specieIndex_];

        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': missing 'file' entry for "
                << BaseModel::modelType() << " model '" << modelType << "'."
                << exit(FatalIOError);
        }

        fileName tableFile = dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': table file \"" << tableFile
                << "\" not found for " << BaseModel::modelType() << "."
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty1D() = default;


    //- Member Functions

    virtual void correct(volScalarField& field) const override
    {
        const objectRegistry& db = this->mesh_;
        tmp<volScalarField> tScalarLookup;

        if (db.template foundObject<volScalarField>(lookupVariable_))
        {
            tScalarLookup = db.template lookupObject<volScalarField>
                (lookupVariable_);
        }
        else if (db.template foundObject<volVectorField>(lookupVariable_))
        {
            tScalarLookup = mag
            (
                db.template lookupObject<volVectorField>(lookupVariable_)
            );
        }
        else
        {
            const word& sName = 
                this->species_.speciesNames()[this->specieIndex_];

            FatalErrorInFunction
                << "Species '" << sName << "': lookupVariable '"
                << lookupVariable_ << "' not found in mesh database!" << nl
                << "Available fields: " << db.names() << exit(FatalError);
        }

        const volScalarField& var = tScalarLookup();

        // 1. Update Internal Field
        field.primitiveFieldRef() = table_.interpolateValues
                                                         (var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            const fvPatchScalarField& pVar = var.boundaryField()[patchi];
            fvPatchScalarField& pField = field.boundaryFieldRef()[patchi];

            forAll(pField, i)
            {
                pField[i] = table_.interpolateValue(pVar[i]);
            }
        }

        field.correctBoundaryConditions();
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty2D Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class TabulatedProperty2D
:
    public BaseModel
{
    //- The 2D interpolation table
    interpolation2DTable<scalar> table_;

    //- The names of the two fields used for lookup (e.g., {"E", "T_e"})
    wordList lookupVariables_;

public:

    //- Constructor
    TabulatedProperty2D
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    )
    :
        BaseModel(modelType, dict, mesh, species, specieIndex),
        table_(dict),
        lookupVariables_(dict.get<wordList>("lookupVariables"))
    {
        const word& sName = this->species_.speciesNames()[this->specieIndex_];

        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': missing 'file' for "
                << BaseModel::modelType() << " model '" << modelType << "'."
                << exit(FatalIOError);
        }

        fileName tableFile = dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': file \"" << tableFile
                << "\" not found for " << BaseModel::modelType() << "."
                << exit(FatalIOError);
        }

        if (lookupVariables_.size() != 2)
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " tabulated2D requires exactly 2 lookupVariables. "
                << "Found: " << lookupVariables_.size() << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty2D() = default;


    //- Member Functions

    virtual void correct(volScalarField& field) const override
    {
        const objectRegistry& db = this->mesh_;

        auto getScalarField = [&](const word& name) -> tmp<volScalarField>
        {
            if (db.template foundObject<volScalarField>(name))
            {
                return db.template lookupObject<volScalarField>(name);
            }
            else if (db.template foundObject<volVectorField>(name))
            {
                return mag(db.template lookupObject<volVectorField>(name));
            }
            const word& sName = this->species_.speciesNames()
                                                        [this->specieIndex_];
            
            FatalErrorInFunction
                << "Species '" << sName << "': lookupVariable '" << name
                << "' not found in mesh database!" << nl
                << "Available fields: " << db.names() << exit(FatalError); 

            return tmp<volScalarField>(nullptr);
        };

        tmp<volScalarField> tX = getScalarField(lookupVariables_[0]);
        tmp<volScalarField> tY = getScalarField(lookupVariables_[1]);

        const volScalarField& fieldX = tX();
        const volScalarField& fieldY = tY();

        forAll(field, cellI)
        {
            field[cellI] = table_(fieldX[cellI], fieldY[cellI]);
        }

        forAll(field.boundaryField(), patchi)
        {
            const fvPatchScalarField& pX = fieldX.boundaryField()[patchi];
            const fvPatchScalarField& pY = fieldY.boundaryField()[patchi];
            fvPatchScalarField& pF = field.boundaryFieldRef()[patchi];

            forAll(pF, i)
            {
                pF[i] = table_(pX[i], pY[i]);
            }
        }

        field.correctBoundaryConditions();
    }
};

/*---------------------------------------------------------------------------*\
                      Class PowerLawProperty Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class PowerLawProperty
:
    public BaseModel
{
    //- Amplitude (A)
    scalar amplitude_;

    //- Exponent (b)
    scalar exponent_;

    //- Lookup field name (e.g., "E")
    word lookupVariable_;

public:

    //- Constructor
    PowerLawProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    )
    :
        BaseModel(modelType, dict, mesh, species, specieIndex),
        amplitude_(0.0),
        exponent_(0.0),
        lookupVariable_("none")
    {
        const word& sName = this->species_.speciesNames()[this->specieIndex_];
        const word ampK = "amplitude";
        const word expK = "exponent";
        const word lookupK = "lookupVariable";

        if (!dict.found(lookupK))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " requires entry '" << lookupK << "'." << exit(FatalIOError);
        }
        dict.readEntry(lookupK, lookupVariable_);

        if (!dict.found(ampK))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " requires entry '" << ampK << "'." << exit(FatalIOError);
        }
        dict.readEntry(ampK, amplitude_);

        if (!dict.found(expK))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " requires entry '" << expK << "'." << exit(FatalIOError);
        }
        dict.readEntry(expK, exponent_);
    }

    //- Member Functions

    virtual void correct(volScalarField& field) const override
    {
        const objectRegistry& db = this->mesh_;

        tmp<volScalarField> tScalarLookup;

        if (db.template foundObject<volScalarField>(lookupVariable_))
        {
            tScalarLookup = db.template lookupObject<volScalarField>
                                                              (lookupVariable_);
        }
        else if (db.template foundObject<volVectorField>(lookupVariable_))
        {
            tScalarLookup = mag
            (
                db.template lookupObject<volVectorField>(lookupVariable_)
            );
        }
        else
        {
            const word& sName = 
                this->species_.speciesNames()[this->specieIndex_];

            FatalErrorInFunction
                << "Species '" << sName << "': lookupVariable '"
                << lookupVariable_ << "' not found in mesh database!" << nl
                << "Available fields: " << db.names() << exit(FatalError);
        }


        const volScalarField& var = tScalarLookup();

        volScalarField dimensionlessVar = 
            var / dimensionedScalar("unit", var.dimensions(), 1.0);

        field = dimensionedScalar("A", field.dimensions(), amplitude_)
              * pow(max(dimensionlessVar, SMALL), exponent_);

        field.correctBoundaryConditions();
    }
};

/*---------------------------------------------------------------------------*\
                      Class Function1Property Declaration
\*---------------------------------------------------------------------------*/

template<class BaseModel>
class Function1Property
:
    public BaseModel
{
    autoPtr<Function1<scalar>> functionPtr_;

    word lookupVariable_;

public:

    //- Constructor
    Function1Property
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex
    )
    :
        BaseModel(modelType, dict, mesh, species, specieIndex),
        functionPtr_(nullptr),
        lookupVariable_("none")
    {
        const word& sName = this->species_.speciesNames()[this->specieIndex_];
        const word fK = "function";

        if (!dict.found("lookupVariable"))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " requires entry 'lookupVariable'." << exit(FatalIOError);
        }
        dict.readEntry("lookupVariable", lookupVariable_);

        if (!dict.found(fK))
        {
            FatalIOErrorInFunction(dict)
                << "Species '" << sName << "': " << BaseModel::modelType()
                << " requires entry '" << fK << "'." << exit(FatalIOError);
        }

        functionPtr_ = Function1<scalar>::New(fK, dict, &mesh);
    }

    //- Member Functions

    virtual void correct(volScalarField& field) const override
    {
        const objectRegistry& db = this->mesh_;
        tmp<volScalarField> tScalarLookup;

        if (db.template foundObject<volScalarField>(lookupVariable_))
        {
            tScalarLookup = db.template lookupObject<volScalarField>
                                                              (lookupVariable_);
        }
        else if (db.template foundObject<volVectorField>(lookupVariable_))
        {
            tScalarLookup = mag(db.template lookupObject<volVectorField>
                                                             (lookupVariable_));
        }
        else
        {
            const word& sName = this->species_.speciesNames()
                                                           [this->specieIndex_];

            FatalErrorInFunction
                << "Species '" << sName << "': lookupVariable '"
                << lookupVariable_ << "' not found in mesh database!" << nl
                << "Available fields: " << db.names() << exit(FatalError);
        }

        const volScalarField& var = tScalarLookup();

        field.primitiveFieldRef() = functionPtr_->value(var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            field.boundaryFieldRef()[patchi] = 
                functionPtr_->value(var.boundaryField()[patchi]);
        }

        field.correctBoundaryConditions();
    }
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
