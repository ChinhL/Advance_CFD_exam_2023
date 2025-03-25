/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = a10fc6a20e1ce523074db5c3c1829c501cd4836e
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void BoundaryLayerInletVelocity_a10fc6a20e1ce523074db5c3c1829c501cd4836e(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    BoundaryLayerInletVelocityMixedValueFvPatchVectorField
);


const char* const BoundaryLayerInletVelocityMixedValueFvPatchVectorField::SHA1sum =
    "a10fc6a20e1ce523074db5c3c1829c501cd4836e";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
BoundaryLayerInletVelocityMixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e"
            " from patch/DimensionedField\n";
    }
}


BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
BoundaryLayerInletVelocityMixedValueFvPatchVectorField
(
    const BoundaryLayerInletVelocityMixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e"
            " from patch/DimensionedField/mapper\n";
    }
}


BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
BoundaryLayerInletVelocityMixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e"
            " from patch/dictionary\n";
    }
}


BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
BoundaryLayerInletVelocityMixedValueFvPatchVectorField
(
    const BoundaryLayerInletVelocityMixedValueFvPatchVectorField& ptf
)
:
    mixedFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e"
            " as copy\n";
    }
}


BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
BoundaryLayerInletVelocityMixedValueFvPatchVectorField
(
    const BoundaryLayerInletVelocityMixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

BoundaryLayerInletVelocityMixedValueFvPatchVectorField::
~BoundaryLayerInletVelocityMixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void BoundaryLayerInletVelocityMixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs BoundaryLayerInletVelocity sha1: a10fc6a20e1ce523074db5c3c1829c501cd4836e\n";
    }

//{{{ begin code
    #line 36 "/home/chinhle/OpenFOAM/chinhle-8/run/examtask/0/U/boundaryField/inlet"
const fvPatch & patch = this->patch();
        const vectorField& Cf=patch.Cf();
        const scalar n=2;
        const scalar Uo=1.5;
        const scalar h=0.4;
        vector flowDir = vector(1,0,0);
       
        forAll (Cf ,faceI)
        {	
        	scalar y=Cf[faceI].y();
        	this->refValue()[faceI]=  Uo*pow(((4*y/h)*(1-(y/h))),n)*flowDir;
        	this->refGrad()[faceI] = vector(0,0,0);//Zero;
        	this->valueFraction()[faceI] = 1;
        	}
//}}} end code

    this->mixedFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

