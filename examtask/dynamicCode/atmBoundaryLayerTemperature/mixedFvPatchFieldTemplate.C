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
    // SHA1 = 4c7c2286d7e05b8b55e0c0964333081e9fc4a553
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void atmBoundaryLayerTemperature_4c7c2286d7e05b8b55e0c0964333081e9fc4a553(bool load)
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
    fvPatchScalarField,
    atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
);


const char* const atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::SHA1sum =
    "4c7c2286d7e05b8b55e0c0964333081e9fc4a553";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553"
            " from patch/DimensionedField\n";
    }
}


atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
(
    const atmBoundaryLayerTemperatureMixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553"
            " from patch/DimensionedField/mapper\n";
    }
}


atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553"
            " from patch/dictionary\n";
    }
}


atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
(
    const atmBoundaryLayerTemperatureMixedValueFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553"
            " as copy\n";
    }
}


atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
atmBoundaryLayerTemperatureMixedValueFvPatchScalarField
(
    const atmBoundaryLayerTemperatureMixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::
~atmBoundaryLayerTemperatureMixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerTemperatureMixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs atmBoundaryLayerTemperature sha1: 4c7c2286d7e05b8b55e0c0964333081e9fc4a553\n";
    }

//{{{ begin code
    #line 36 "/home/chinhle/OpenFOAM/chinhle-8/run/examtask/0/T/boundaryField/bottom"
const fvPatch & patch = this->patch();
        const vectorField& Cf=patch.Cf();
       
        forAll (Cf ,faceI)
        {
        	scalar x=Cf[faceI].x();
        	if((Foam::sin((constant::mathematical::pi*x)/1))>0)
        	{ 
        	this->refValue()[faceI]=  400;
        	this->refGrad()[faceI] = 0;//Zero;
        	this->valueFraction()[faceI] = 1.0;
        	}
        	else { 
        	this->refValue()[faceI]=  0;
        	this->refGrad()[faceI] = 0;//Zero;
        	this->valueFraction()[faceI] = 1.0;
        	}
        	}
//}}} end code

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

