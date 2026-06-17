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
    Foam::photoionizationModel
    
Description
    Abstract base class for photoionization source term models.

    Owns the photoionization source field Sph [1/(m3 s)], common to all
    concrete models. Derived classes implement correct() to recompute Sph
    given the electron-impact ionization source Si.

SourceFiles
    photoionizationModel.C

\*---------------------------------------------------------------------------*/

#ifndef photoionizationModel_H
#define photoionizationModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class photoionizationModel Declaration
\*---------------------------------------------------------------------------*/

class photoionizationModel
:
    public IOdictionary
{
private:

    //- Disallow copy
    photoionizationModel(const photoionizationModel&) = delete;
    void operator=(const photoionizationModel&) = delete;

protected:

    //- Reference to the primary gas/plasma mesh
    const fvMesh& mesh_;

    //- Photoionization source [1/(m3 s)]
    volScalarField Sph_;

    //- Interval of solving for photoionization
    label solveInterval_;

public:

    //- Runtime type information
    TypeName("photoionizationModel");

    //- Declare the Run-Time Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        photoionizationModel,
        dictionary,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        photoionizationModel(const fvMesh& mesh);

    // Selectors

        //- Selector
        static autoPtr<photoionizationModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~photoionizationModel() = default;

    // Member Functions

        //- Update the photoionization source field
        virtual void correct() = 0;

        //- Access the photoionization source field
        const volScalarField& Sph() const noexcept { return Sph_; }

        //- Access the solve interval
        label solveInterval() const noexcept { return solveInterval_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
