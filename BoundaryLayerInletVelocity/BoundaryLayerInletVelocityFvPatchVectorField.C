/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "BoundaryLayerInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::BoundaryLayerInletVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoundaryLayerInletVelocityFvPatchVectorField::
BoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowDir_(Zero),
    yDir_(Zero),
    n_(2),
    Uo_(1.5),
    h_(0.4)
{
}


Foam::BoundaryLayerInletVelocityFvPatchVectorField::
BoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    flowDir_(dict.lookup<vector>("flowDir")),
    yDir_(dict.lookup<vector>("yDir")),
    n_(dict.lookupOrDefault<scalar>("n", 2)),
    Uo_(dict.lookupOrDefault<scalar>("Uo", 1.5)),
    h_(dict.lookupOrDefault<scalar>("h", 0.4))
{


    fixedValueFvPatchVectorField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
    */
}


Foam::BoundaryLayerInletVelocityFvPatchVectorField::
BoundaryLayerInletVelocityFvPatchVectorField
(
    const BoundaryLayerInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    flowDir_(ptf.flowDir_),
    yDir_(ptf.yDir_),
    n_(ptf.n_),
    Uo_(ptf.Uo_),
    h_(ptf.h_)
{}


Foam::BoundaryLayerInletVelocityFvPatchVectorField::
BoundaryLayerInletVelocityFvPatchVectorField
(
    const BoundaryLayerInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    flowDir_(ptf.flowDir_),
    yDir_(ptf.yDir_),
    n_(ptf.n_),
    Uo_(ptf.Uo_),
    h_(ptf.h_)
{}


Foam::BoundaryLayerInletVelocityFvPatchVectorField::
BoundaryLayerInletVelocityFvPatchVectorField
(
    const BoundaryLayerInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    flowDir_(ptf.flowDir_),
    yDir_(ptf.yDir_),
    n_(ptf.n_),
    Uo_(ptf.Uo_),
    h_(ptf.h_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BoundaryLayerInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    //m(fieldData_, fieldData_);
}


void Foam::BoundaryLayerInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const BoundaryLayerInletVelocityFvPatchVectorField& tiptf =
        refCast<const BoundaryLayerInletVelocityFvPatchVectorField>(ptf);

    //fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::BoundaryLayerInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    /*fixedValueFvPatchVectorField::operator==
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );*/
    const fvPatch& patch= this->patch();
    const vectorField& Cf= patch.Cf();
    const scalarField YY =(Cf & yDir_);
    scalarField U=Uo_*pow(((4*YY/h_)*(1-(YY/h_))),n_);

    fixedValueFvPatchVectorField::operator==(flowDir_*U);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::BoundaryLayerInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "flowDir", flowDir_);
    writeEntry(os, "yDir", yDir_);
    writeEntry(os, "Uo", Uo_);
    writeEntry(os, "n", n_);
    writeEntry(os, "h", h_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        BoundaryLayerInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
