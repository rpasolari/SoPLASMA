/*---------------------------------------------------------------------------*\
  File: updateSurfaceCharge.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Updates the global surface charge density field surfCharge by integrating 
    the net wall current density over the current time-step. 

    The update iterates over all boundary patches and all species. It 
    accumulates contributions only from boundaries where the 
    'solidSurfaceFlux' condition is active and the 'enableSurfaceCharging' 
    switch is enabled by the user.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info << "Updating surface charge!" << endl;

const fvBoundaryMesh& boundary = gasMesh().boundary();

forAll(boundary, patchi)
{
    const fvPatch& p = boundary[patchi];

    bool isDielectricInterface = false;
    word neighborName = "none";

    // Chech if the wall is type mapped
    if (isA<mappedPatchBase>(p.patch()))
    {
        const mappedPatchBase& mpp =
            refCast<const mappedPatchBase>(p.patch());
        neighborName = mpp.sampleRegion();

        // Check if region is one of the dielectric regions
        forAll(dielectricNames, d)
        {
            if (neighborName == dielectricNames[d])
            {
                isDielectricInterface = true;
                break;
            }
        }
    }

    // Proceed only if it is a dielectric wall
    if (!isDielectricInterface)
    {
        continue;
    }

    scalarField& sigmaPatch = surfCharge.boundaryFieldRef()[patchi];

    // Iterate over all species to check their boundary conditions
    for (label i = 0; i < species.nSpecies(); ++i)
    {
        const volScalarField& speciesField = species.numberDensity(i);
        const fvPatchField<scalar>& pField =
            speciesField.boundaryField()[patchi];

        // Check if it is a solidSurfaceFlux (or its derived b.c.) type
        if (isA<ddSolidSurfaceFluxFvPatchScalarField>(pField))
        {
            const auto& bc = 
                refCast<const ddSolidSurfaceFluxFvPatchScalarField>(pField);
            
            // Check the switch for adding or not to surfCharge
            bool writeOn = bc.enableSurfaceCharging();
            
            if (writeOn)
            {
                const scalar dt = runTime.deltaT().value();
                const scalarField& magSf =
                    gasMesh().magSf().boundaryField()[patchi];

                // Get flux
                const scalarField& phi_p = 
                    transport.particleFlux(i).boundaryField()[patchi];

                // Get charge
                scalar q = species.speciesCharge(i).value();

                // Add to surfCharge
                // sigmaPatch += (phi_p * q * dt) / magSf;
                sigmaPatch = sigma0.boundaryField()[patchi] + (phi_p * q * dt) / magSf;
            }
        }
    }
}
