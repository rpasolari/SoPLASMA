/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
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
#include "interpolationTable.H"
#include "interpolation2DTable.H"

namespace Foam
{

// Used by all non-constant templates

namespace detail
{

template<class Mesh>
tmp<volScalarField> getScalarField
(
    const Mesh& mesh,
    const word& name,
    const word& propertyName
)
{
    if (mesh.template foundObject<volScalarField>(name))
    {
        return mesh.template lookupObject<volScalarField>(name);
    }

    if (mesh.template foundObject<volVectorField>(name))
    {
        return mag(mesh.template lookupObject<volVectorField>(name));
    }

    const wordList availableScalars =
        mesh.template names<volScalarField>();
    const wordList availableVectors =
        mesh.template names<volVectorField>();

    FatalErrorInFunction
        << "plasmaPropertyEvaluator '" << propertyName
        << "': lookupVariable '" << name
        << "' not found in mesh registry." << nl
        << "Available volScalarFields: " << availableScalars << nl
        << "Available volVectorFields: " << availableVectors << nl
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}

} // End namespace detail

/*---------------------------------------------------------------------------*\
                      Class ConstantProperty Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class ConstantProperty
:
    public Base
{
    scalar val_;

public:
    ConstantProperty
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        val_(0.0)
    {
        if (!dict.readIfPresent("value", val_))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': constant model requires entry 'value'." << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~ConstantProperty() = default;

    //- Member Functions

    bool isUniform() const override { return true; }
    bool isConstant() const override { return true; }

    dimensionedScalar evaluateScalar() const override
    {
        return dimensionedScalar(this->propertyName(), this->dims_, val_);
    }

    tmp<volScalarField> evaluate() const override
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, val_)
        );
    }

    void correct(volScalarField& field) const override
    {
        field.primitiveFieldRef() = val_;

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                field.boundaryFieldRef()[patchi] = val_;
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch(scalarField& pField, const label) const override
    {
        pField = val_;
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty1D Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class TabulatedProperty1D
:
    public Base
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
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        table_(dict),
        lookupVariable_(dict.get<word>("lookupVariable"))
    {
        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated1D requires entry 'file'." << nl
                << exit(FatalIOError);
        }

        const fileName tableFile =
            dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': table file not found:" << nl
                << "    " << tableFile << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty1D() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on tabulated1D evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            table_.interpolateValues(var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                evaluatePatch(field.boundaryFieldRef()[patchi], patchi);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        forAll(pField, i)
        {
            pField[i] = table_.interpolateValue(pVar[i]);
        }
    }
};

/*---------------------------------------------------------------------------*\
                      Class TabulatedProperty2D Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class TabulatedProperty2D
:
    public Base
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
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        table_(dict),
        lookupVariables_(dict.get<wordList>("lookupVariables"))
    {
        if (lookupVariables_.size() != 2)
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated2D requires exactly 2 lookupVariables." << nl
                << exit(FatalIOError);
        }

        if (!dict.found("file"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': tabulated2D requires entry 'file'." << nl
                << exit(FatalIOError);
        }

        const fileName tableFile =
            dict.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': table file not found:" << nl
                << "    " << tableFile << nl
                << exit(FatalIOError);
        }
    }

    //- Destructor
    virtual ~TabulatedProperty2D() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on tabulated2D evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tX = detail::getScalarField
        (
            this->mesh_, lookupVariables_[0], this->propertyName()
        );
        tmp<volScalarField> tY = detail::getScalarField
        (
            this->mesh_, lookupVariables_[1], this->propertyName()
        );

        const volScalarField& X = tX();
        const volScalarField& Y = tY();

        forAll(field.primitiveField(), cellI)
        {
            field.primitiveFieldRef()[cellI] = table_(X[cellI], Y[cellI]);
        }

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                evaluatePatch(field.boundaryFieldRef()[patchi], patchi);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tX = detail::getScalarField
        (
            this->mesh_, lookupVariables_[0], this->propertyName()
        );
        tmp<volScalarField> tY = detail::getScalarField
        (
            this->mesh_, lookupVariables_[1], this->propertyName()
        );

        const scalarField& pX = tX().boundaryField()[patchi];
        const scalarField& pY = tY().boundaryField()[patchi];

        forAll(pField, i)
        {
            pField[i] = table_(pX[i], pY[i]);
        }
    }
};

/*---------------------------------------------------------------------------*\
                      Class PowerLawProperty Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class PowerLawProperty
:
    public Base
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
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        amplitude_(0.0),
        exponent_(0.0),
        lookupVariable_("none")
    {
        if (!dict.found("amplitude"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'amplitude'." << nl
                << exit(FatalIOError);
        }
        amplitude_ = dict.get<scalar>("amplitude");

        if (!dict.found("exponent"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'exponent'." << nl
                << exit(FatalIOError);
        }
        exponent_ = dict.get<scalar>("exponent");

        if (!dict.found("lookupVariable"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': powerLaw requires entry 'lookupVariable'." << nl
                << "Provide the name of the field to use as the independent "
                << "variable (e.g. 'Emag', 'reducedE', 'Te')." << nl
                << exit(FatalIOError);
        }
        lookupVariable_ = dict.get<word>("lookupVariable");
    }

    //- Destructor
    virtual ~PowerLawProperty() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on powerLaw evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            amplitude_ * pow(max(var.primitiveField(), SMALL), exponent_);

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                evaluatePatch(field.boundaryFieldRef()[patchi], patchi);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        pField = amplitude_ * pow(max(pVar, SMALL), exponent_);
    }
};

/*---------------------------------------------------------------------------*\
                      Class Function1Property Declaration
\*---------------------------------------------------------------------------*/

template<class Base>
class Function1Property
:
    public Base
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
        const word& propertyName,
        const dimensionSet& dims
    )
    :
        Base(modelType, dict, mesh, propertyName, dims),
        functionPtr_(nullptr),
        lookupVariable_("none")
    {
        if (!dict.found("lookupVariable"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': function1 requires entry 'lookupVariable'." << nl
                << "Provide the name of the field to use as the independent "
                << "variable (e.g. 'Emag', 'reducedE', 'Te')." << nl
                << exit(FatalIOError);
        }
        lookupVariable_ = dict.get<word>("lookupVariable");

        if (!dict.found("function"))
        {
            FatalIOErrorInFunction(dict)
                << "plasmaPropertyEvaluator '" << propertyName
                << "': function1 requires entry 'function'." << nl
                << "Provide a valid OpenFOAM Function1<scalar> entry." << nl
                << exit(FatalIOError);
        }
        functionPtr_ = Function1<scalar>::New("function", dict, &mesh);
    }

    //- Destructor
    virtual ~Function1Property() = default;

    //- Member Functions

    bool isUniform() const override { return false; }
    bool isConstant() const override { return false; }

    dimensionedScalar evaluateScalar() const override
    {
        FatalErrorInFunction
            << "evaluateScalar() called on function1 evaluator '"
            << this->propertyName() << "' which is not uniform." << nl
            << exit(FatalError);
        return dimensionedScalar(this->propertyName(), this->dims_, 0.0);
    }

    tmp<volScalarField> evaluate() const override
    {
        auto tField = tmp<volScalarField>::New
        (
            IOobject
            (
                this->propertyName(),
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensionedScalar(this->propertyName(), this->dims_, 0.0)
        );
        correct(tField.ref());
        return tField;
    }

    void correct(volScalarField& field) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const volScalarField& var = tVar();

        field.primitiveFieldRef() =
            functionPtr_->value(var.primitiveField());

        forAll(field.boundaryField(), patchi)
        {
            if (!field.boundaryField()[patchi].coupled())
            {
                evaluatePatch(field.boundaryFieldRef()[patchi], patchi);
            }
        }

        field.correctBoundaryConditions();
    }

    void evaluatePatch
    (
        scalarField& pField,
        const label patchi
    ) const override
    {
        tmp<volScalarField> tVar = detail::getScalarField
        (
            this->mesh_, lookupVariable_, this->propertyName()
        );
        const scalarField& pVar = tVar().boundaryField()[patchi];

        pField = functionPtr_->value(pVar);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
