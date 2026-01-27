/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
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
    Foam::decompositionConstraints::preserveMappedInterface

Description
    Constraint to synchronize processor assignment across mapped interfaces.

    Ensures that cells adjacent to mapped patches remain on the same processor
    across different regions. This is critical for solvers where 1-to-1 face
    matching is required or where neighbor data must be accessed without
    inter-processor communication overhead.

Usage:
    \verbatim
    constraints
    {
        syncInterface
        {
            type    preserveMappedInterface;
            patches (interfacePatchName);
        }
    }
    \endverbatim

SourceFiles
    preserveMappedInterfaceConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef preserveFaceZonesConstraint_H
#define preserveFaceZonesConstraint_H

#include "decompositionConstraint.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{

/*---------------------------------------------------------------------------*\
                      Class preserveMappedInterface Declaration
\*---------------------------------------------------------------------------*/

class preserveMappedInterface
:
    public decompositionConstraint
{
    // Private data

        //- List of patches to keep together
        wordRes patches_;


public:

    //- Runtime type information
    TypeName("preserveMappedInterface");


    // Constructors

        //- Construct with constraint dictionary
        explicit preserveMappedInterface(const dictionary& dict);

        //- Construct from components
        explicit preserveMappedInterface(const UList<wordRe>& zones);


    //- Destructor
    virtual ~preserveMappedInterface() = default;


    // Member Functions

        //- Add my constraints to list of constraints
        virtual void add
        (
            const polyMesh& mesh,
            boolList& blockedFace,
            PtrList<labelList>& specifiedProcessorFaces,
            labelList& specifiedProcessor,
            List<labelPair>& explicitConnections
        ) const;

        //- Apply any additional post-decomposition constraints
        virtual void apply
        (
            const polyMesh& mesh,
            const boolList& blockedFace,
            const PtrList<labelList>& specifiedProcessorFaces,
            const labelList& specifiedProcessor,
            const List<labelPair>& explicitConnections,
            labelList& decomposition
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionConstraints
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
