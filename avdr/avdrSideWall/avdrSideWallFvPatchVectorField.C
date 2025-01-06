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

\*---------------------------------------------------------------------------*/

#include "avdrSideWallFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::avdrSideWallFvPatchVectorField::
avdrSideWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    avdr_(1.0)
{
}


Foam::avdrSideWallFvPatchVectorField::
avdrSideWallFvPatchVectorField
(
    const avdrSideWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    avdr_(ptf.avdr_)
{}


Foam::avdrSideWallFvPatchVectorField::
avdrSideWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    avdr_(readScalar(dict.lookup("avdr")))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::avdrSideWallFvPatchVectorField::
avdrSideWallFvPatchVectorField
(
    const avdrSideWallFvPatchVectorField& sfspvf
)
  :
  fixedValueFvPatchVectorField(sfspvf),
  avdr_(sfspvf.avdr_)
{}


Foam::avdrSideWallFvPatchVectorField::
avdrSideWallFvPatchVectorField
(
    const avdrSideWallFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
  fixedValueFvPatchVectorField(sfspvf, iF),
  avdr_(sfspvf.avdr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::avdrSideWallFvPatchVectorField::updateCoeffs()
{
  if (updated())
  {
    return;
  }

  const surfaceScalarField& phi =
   db().lookupObject<surfaceScalarField>("phi");

  const volVectorField& U = 
   db().lookupObject<volVectorField>("U");

  const volScalarField& rho = 
    db().lookupObject<volScalarField>("rho");

  const fvMesh & mesh = U.mesh();

  label  inletId  = mesh.boundaryMesh().findPatchID("inlet" );
  label outletId  = mesh.boundaryMesh().findPatchID("outlet");

  const scalar m_dotInlet  = -gSum(phi.boundaryField()[inletId] );
  const scalar m_dotOutlet =  gSum(phi.boundaryField()[outletId]);  

  Info << "AVDR = " << m_dotOutlet/m_dotInlet << endl;  
  const scalar avdr = avdr_;

  const scalar m_dotAddOnePatch = (avdr-1)*m_dotInlet/2;

  volScalarField distribution
  (
    IOobject
  (
    "distribution",
    "system",
    mesh,
    IOobject::MUST_READ, //READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
    mesh
  );

  const labelUList& owner = mesh.owner();
      
  scalar aAVDR = 0.;
  label patchName = patch().index();
  vectorField n = patch().nf();

  for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); ++faceI)
  {
    label onThisPatch = mesh.boundaryMesh().whichPatch(faceI);
    if(patchName == onThisPatch)
    {
      label own = owner[faceI];
      aAVDR = aAVDR + mesh.magSf().boundaryField()[onThisPatch][faceI-mesh.boundary()[onThisPatch].start()]*distribution[own];
    }
  }

  const scalar rhoV = m_dotAddOnePatch/aAVDR;

  fvPatchField<vector> patchU = U.boundaryField()[patchName];

  for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); ++faceI)
  {
    label onThisPatch = mesh.boundaryMesh().whichPatch(faceI);
    if(patchName == onThisPatch)
    {
      label own = owner[faceI];
      vector Ui
      (
        U[own][0], 
        U[own][1], 
        -sign(n[faceI-mesh.boundary()[onThisPatch].start()][2])*rhoV/rho[own]*distribution[own]
      );
      patchU[faceI-mesh.boundary()[onThisPatch].start()] = Ui;
    }
  }
  operator==(patchU);
  fixedValueFvPatchVectorField::updateCoeffs(); 
}

void Foam::avdrSideWallFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("avdr", avdr_);

    this->writeEntry("value", os);

}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::avdrSideWallFvPatchVectorField::operator=
(
 const fvPatchField<vector>& patchU
 )
{
  fvPatchField<vector>::operator=(patchU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        avdrSideWallFvPatchVectorField
    );
}

// ************************************************************************* //