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
    Foam::noPhotoionization
    
Description
    Null photoionization model. The source term remains zero throughout the
    simulation.

    Selected in constant/photoionizationProperties via

        photoionizationModel  none;

SourceFiles
    noPhotoionization.C

\*---------------------------------------------------------------------------*/

#ifndef noPhotoionization_H
#define noPhotoionization_H

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class noPhotoionization Declaration
\*---------------------------------------------------------------------------*/

class noPhotoionization
:
    public photoionizationModel
{
private:

    //- Disallow copy
    noPhotoionization(const noPhotoionization&) = delete;
    void operator=(const noPhotoionization&) = delete;

public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from mesh
        noPhotoionization(const fvMesh& mesh);

    //- Destructor
    virtual ~noPhotoionization() = default;

    // Member Functions

        //- Return the photoionization source term (always zero)
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
