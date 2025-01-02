/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::avdrSideWall

Description
    - should be used with zero gradient for all other variables

Example

    front
    {
        type            avdrSideWall;
        avdr            1.1;
        value           $internalField;
    } 

\*---------------------------------------------------------------------------*/

#ifndef avdrSideWallFvPatchVectorFields_H
#define avdrSideWallFvPatchVectorFields_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class avdrSideWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class avdrSideWallFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

    // TODO: add inlet and outlet as parameters with default

    //- required AVDR
    scalar avdr_;

public:

    //- Runtime type information
    TypeName("avdrSideWall");


    // Constructors

        //- Construct from patch and internal field
        avdrSideWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        avdrSideWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given avdrSideWallFvPatchVectorField
        //  onto a new patch
        avdrSideWallFvPatchVectorField
        (
            const avdrSideWallFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        avdrSideWallFvPatchVectorField
        (
            const avdrSideWallFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new avdrSideWallFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        avdrSideWallFvPatchVectorField
        (
            const avdrSideWallFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new avdrSideWallFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const scalar& avdr() 
            {
                return avdr_;
            }

    virtual void updateCoeffs();
    
    
    //- Write
    virtual void write(Ostream&) const;
    
    virtual void operator=(const fvPatchField<vector>& pvf);
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //