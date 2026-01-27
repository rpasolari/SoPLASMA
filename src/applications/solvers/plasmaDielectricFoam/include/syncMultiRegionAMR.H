
/*---------------------------------------------------------------------------*\
  File: syncMultiRegionAMR.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Synchronizes refinement signals across coupled regions and manages 
    topological updates. Includes logic to reset LDU matrix assembly and 
    clear mapping caches to maintain implicit coupling stability.

    The multi-region refinement logic implemented in this toolkit is based on 
    the methodology developed in the blastAMR project:
    https://github.com/STFS-TUDa/blastAMR

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if 
(
    enableAMR &&
    totalRegions > 1 &&
    (
        (runTime.timeIndex() % refineInterval == 0) 
     || (runTime.timeIndex() % unrefineInterval == 0)
     || (runTime.timeIndex() % balanceInterval == 0)
    )
)
{
    // Read the refinement field of the gas region, and build the refineField
    if
    (
        !restrictAMRToPlasma
        && conformalAMRSync 
        && refineFieldGasPtr.valid()
    )
    {
        gasAMRMesh.amr().error().update(); 

        volScalarField& refineFieldGas = refineFieldGasPtr();
        
        refineFieldGas.primitiveFieldRef() = gasAMRMesh.amr().error().error();
        refineFieldGas.correctBoundaryConditions();

        forAll(dielectricRegions, i)
        {
            if (touchesGas[i] && refineFieldDielectricList.set(i))
            {
                volScalarField& refineFieldDiel = refineFieldDielectricList[i];
                refineFieldDiel.primitiveFieldRef() = 0.0;
                refineFieldDiel.correctBoundaryConditions();
            }
        }
    }
}

if 
(
    enableAMR &&
    (
        (runTime.timeIndex() % refineInterval == 0) 
     || (runTime.timeIndex() % unrefineInterval == 0)
     || (runTime.timeIndex() % balanceInterval == 0)
    )
)
{
    bool gasChanged = gasAMRMesh.update();

    if (totalRegions > 1)
    {
        bool dielectricChanged = false;
        forAll(dielectricRegions, i)
        {
            adaptiveFvMesh& dielAMRMesh = 
                refCast<adaptiveFvMesh>(dielectricRegions[i]);

            if (dielAMRMesh.update())
            {
                dielectricChanged = true;
            }
        }
    

        // If any mesh updated cleanup
        if (gasChanged || dielectricChanged)
        {

            auto clearMappedCache = [](const fvMesh& mesh)
            {
                forAll(mesh.boundaryMesh(), patchi)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchi];
                    const mappedPatchBase* mpPtr = 
                        dynamic_cast<const mappedPatchBase*>(&pp);

                    if (mpPtr)
                    {
                        const_cast<mappedPatchBase*>(mpPtr)->clearOut();
                    }
                }
            };

            clearMappedCache(gasMesh());

            forAll(dielectricRegions, i)
            {
                clearMappedCache(dielectricRegions[i]);
            }

            // If coupled (implicit coupling) reset matrix assembly
            if (coupled)
            {
                gasMesh().thisDb().checkOut("lduAssembly4");

                fvMatrixAssemblyPtr.reset(nullptr);
                fvMatrixAssemblyPtr.reset
                (
                    new fvMatrix<scalar>
                    (
                        ePotential,
                        dimensionSet(0, 0, 1, 0, 0, 1, 0)
                    )
                );
            }            
        }
    }
}
