/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::coupledElectricPotentialFvPatchScalarField

Description
    Based on compressible::turbulentTemperatureRadCoupledMixedFvPatchScalarField
    in OpenFOAM v2412, adapted for electrostatic plasma potential coupling.

    Coupled boundary condition for electric potential between adjacent regions,
    enforcing potential continuity and optional surface charge effects at
    plasma–dielectric or plasma–solid interfaces.

    The condition ensures:
        - Potential continuity across the interface
        - Normal electric flux jump for surface charge density sigmaS

    The coupling uses an implicit formulation to ensure consistent field
    coupling between regions.

Usage
    Example boundary condition entry:
    \verbatim
    <patchName>
    {
    type            coupledElectricPotential;

    // Name of electric potential field in the neighbor region
    phiNbr          ePotential;

    // Only used for restart / initial reference in mixed BC
    value           uniform 0;

    // If true, patches are coupled, and the equation is solved in a 
    // monolithic way
    useImplicit     true;

    // Optional: write interface diagnostics to log
    verbose         true;

    // Optional: write report every N time steps (requires verbose=true)
    logInterval     1;

    // Surface charge density [C/m^2] on this side of the interface.
    // Must correspond to a volScalarField named 'surfCharge' in this 
    // region’s 0/ directory. Set to "none" to disable.
    surfCharge      surfCharge;

    // Surface charge density [C/m^2] on the neighbor side. Must exist as
    // a volScalarField in the neighbor region with the same patch name.
    // Set to "none" if the neighbor region has zero surface charge.
    surfChargeNbr   none;
    }
    \endverbatim

Notes
    - Works only on mapped, AMI, or cyclic-type interface patches.
    - In parallel runs, AMI patch-pairs must reside on the same processor.
    - For cyclic and mapped patches, the number of faces per processor must 
       match (requires a constraint-aware decomposition method).    

SourceFiles
    coupledElectricPotentialFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_coupledElectricPotentialFvPatchScalarField_H
#define Foam_coupledElectricPotentialFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarList.H"
#include "mappedPatchFieldBase.H"
#include "writeFile.H"
#include "scalarField.H"
#include "Enum.H"
#include "fvPatch.H"
#include "PatchFunction1.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "IOdictionary.H" 

#include "electromagneticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class coupledElectricPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledElectricPotentialFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public mappedPatchFieldBase<scalar>,
    public functionObjects::writeFile
{
    // Private Data

        //- Name of potential field on the neighbour region
        const word phiNbrName_;

        //- Name of the surface charge in the neighbour region
        const word surfChargeNbrName_;

        //- Name of the surface charge in this region
        const word surfChargeName_;
    
        //- Log-file output frequency [s]
        scalar logInterval_;

        //- Index of last write
        label executionIndex_;

        //- Flag to enable verbose output with debug information
        const bool verbose_;

        //- Prefix for auxiliary field names (in verbose mode)
        const word prefix_;


    // Private Member Functions

        //- Return delta phi between regions
        tmp<scalarField> deltaPhi() const;

        //- Return the sum of deltaCoeff*epsilon from nbr and local
        tmp<scalarField> beta() const;

        //- Return requested field from the object registry
        //- or create+register the field to the object registry
        volScalarField& getOrCreateField(const word& fieldName) const;

        //- Store electric related fields to facilitate postprocessing
        void storeECFields
        (
            const word& prefix,
            const scalarField& etc,
            const scalarField& etcPatch
        ) const;

        //- Flag to allow writing to file
        bool writeFile();

        // functionObjects::writeFile

            //- Write file header information
            void writeFileHeader(Ostream& os);

public:

    //- Runtime type information
    TypeName("coupledElectricPotential");


    // Constructors

        //- Construct from patch and internal field
        coupledElectricPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledElectricPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        coupledElectricPotentialFvPatchScalarField
        (
            const coupledElectricPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        coupledElectricPotentialFvPatchScalarField
        (
            const coupledElectricPotentialFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        coupledElectricPotentialFvPatchScalarField
        (
            const coupledElectricPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map given fvPatchField onto this fvPatchField
        virtual void rmap(const fvPatchField<scalar>&, const labelList&);


    // Evaluation

        //- Compute absolute permittivity ε = ε0 * εr 
        scalar epsilon() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Manipulate matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& m,
            const label iMatrix,
            const direction cmpt
        );

        virtual void initEvaluate
        (
            const UPstream::commsTypes commsType =
                UPstream::commsTypes::blocking
        ) override;

        //- Final evaluation of patch values after matrix solve
        //- Ensures mixed-type update is applied using the computed coefficients
        virtual void evaluate
        (
            const Pstream::commsTypes commsType =
                Pstream::commsTypes::buffered
        );

    // I-O

        //- Write dictionary entries
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
