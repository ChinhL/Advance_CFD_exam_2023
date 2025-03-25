/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "temperatureLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(temperatureLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        temperatureLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::temperatureLaw::calcNu() const
{
    const volScalarField& T_=U_.mesh().lookupObject<volScalarField>("T");
    return min
    (
        nuMax_,
        max
        (
            nuMin_,
            nuRef_-(T_-TRef_)*nuDeg_
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::temperatureLaw::temperatureLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    temperatureLawCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    TRef_("TRef", dimensionSet (0 ,0, 0 ,1 ,0 ,0 ,0), temperatureLawCoeffs_),
    
    nuMin_("nuMin", dimViscosity, temperatureLawCoeffs_),
    nuMax_("nuMax", dimViscosity, temperatureLawCoeffs_),
    nuRef_("nuRef", dimViscosity, temperatureLawCoeffs_),
    nuDeg_("nuDeg", dimensionSet (0 ,2, -1 ,-1 ,0 ,0 ,0), temperatureLawCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::temperatureLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    temperatureLawCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    temperatureLawCoeffs_.lookup("TRef") >> TRef_;
    temperatureLawCoeffs_.lookup("nuMin") >> nuMin_;
    temperatureLawCoeffs_.lookup("nuMax") >> nuMax_;
    temperatureLawCoeffs_.lookup("nuRef") >> nuRef_;
    temperatureLawCoeffs_.lookup("nuDeg") >> nuDeg_;


    return true;
}


// ************************************************************************* //
