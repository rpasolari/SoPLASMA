/*---------------------------------------------------------------------------*\
  File: preserveMappedInterfaceConstraint.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::decompositionConstraints::preserveMappedInterface

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "preserveMappedInterfaceConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "mappedPolyPatch.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preserveMappedInterface);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preserveMappedInterface,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preserveMappedInterface::preserveMappedInterface
(
    const dictionary& dict
)
:
    decompositionConstraint(dict, typeName),
    patches_(coeffDict_.get<wordRes>("patches"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : locking processor assignment for patches " 
            << patches_ << endl;
    }
}


Foam::decompositionConstraints::preserveMappedInterface::preserveMappedInterface
(
    const UList<wordRe>& patches
)
:
    decompositionConstraint(dictionary(), typeName),
    patches_(patches)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : locking processor assignment for patches " 
            << patches_ << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preserveMappedInterface::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelList patchIDs(patches_.matching(pbm.names()));

    label nUnblocked = 0;

    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, i)
        {
            label meshFacei = pp.start() + i;
            if (blockedFace[meshFacei])
            {
                blockedFace[meshFacei] = false;
                ++nUnblocked;
            }
        }
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::preserveMappedInterface::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const labelList patchIDs(patches_.matching(pbm.names()));
    
    if (patchIDs.empty())
    {
        Info<< " No matching patches found in mesh " << mesh.name() 
            << ". Check your decomposeParDict names!" << endl;
    }

    label nChanged = 0;
    do
    {
        nChanged = 0;
        for (const label patchi : patchIDs)
        {
            // Check for both generic mapped patches and mappedWall patches
            if 
            (
                isType<mappedPolyPatch>(pbm[patchi]) 
             || isType<mappedWallPolyPatch>(pbm[patchi])
            )
            {
                const polyPatch& pp = pbm[patchi];
                const labelList& faceCells = pp.faceCells();
                
                for (const label celli : faceCells)
                {
                    if (decomposition[celli] != Pstream::myProcNo())
                    {
                        decomposition[celli] = Pstream::myProcNo();
                        ++nChanged;
                    }
                }
            }
        }

    } while (returnReduceOr(nChanged));

}


// ************************************************************************* //
