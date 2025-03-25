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

#include "atmBoundaryLayerTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::atmBoundaryLayerTemperatureFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmBoundaryLayerTemperatureFvPatchScalarField::
atmBoundaryLayerTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    scalarData_(0.0),
    xDir_(Zero),
    fieldData_(p.size(), Zero),
    timeVsData_(),
    wordData_("wordDefault"),
    labelData_(-1),
    boolData_(false)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::atmBoundaryLayerTemperatureFvPatchScalarField::
atmBoundaryLayerTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    scalarData_(dict.lookup<scalar>("scalarData")),
    xDir_(dict.lookup<vector>("xDir")),
    fieldData_("fieldData", dict, p.size()),
    timeVsData_(Function1<scalar>::New("timeVsData", dict)),
    wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
    labelData_(-1),
    boolData_(false)
{
    refGrad() = Zero;
    valueFraction() = 0.0;

    refValue() = scalarField("value", dict, p.size());
    fvPatchScalarField::operator=(refValue());


    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    refValue() = *this;
    */
}


Foam::atmBoundaryLayerTemperatureFvPatchScalarField::
atmBoundaryLayerTemperatureFvPatchScalarField
(
    const atmBoundaryLayerTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    xDir_(ptf.xDir_),
    fieldData_(mapper(ptf.fieldData_)),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


Foam::atmBoundaryLayerTemperatureFvPatchScalarField::
atmBoundaryLayerTemperatureFvPatchScalarField
(
    const atmBoundaryLayerTemperatureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    scalarData_(ptf.scalarData_),
    xDir_(ptf.xDir_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


Foam::atmBoundaryLayerTemperatureFvPatchScalarField::
atmBoundaryLayerTemperatureFvPatchScalarField
(
    const atmBoundaryLayerTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    scalarData_(ptf.scalarData_),
    xDir_(ptf.xDir_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmBoundaryLayerTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(fieldData_, fieldData_);
}


void Foam::atmBoundaryLayerTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const atmBoundaryLayerTemperatureFvPatchScalarField& tiptf =
        refCast<const atmBoundaryLayerTemperatureFvPatchScalarField>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::atmBoundaryLayerTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    /*mixedFvPatchScalarField::refValue() =
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "phi"
        );
    valueFraction() = 1.0 - pos0(phip);*/
    const fvPatch& patch= this->patch();
    const vectorField& Cf= patch.Cf();
    const scalarField XX =(Cf & xDir_);
    //if((Foam::sin((constant::mathematical::pi*XX)/0.1))>0.0) {const double T=400;}
    const double T=400;
    mixedFvPatchScalarField::refValue()=T;	
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::atmBoundaryLayerTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "scalarData", scalarData_);
    writeEntry(os, "xDir", xDir_);
    writeEntry(os, "fieldData", fieldData_);
    writeEntry(os, timeVsData_());
    writeEntry(os, "wordData", wordData_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        atmBoundaryLayerTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
