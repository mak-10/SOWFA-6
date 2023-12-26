/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "windFarmCanopy.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(windFarmCanopy, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::windFarmCanopy::read()
{
    dictionary subDict = dict_.subOrEmptyDict("windFarmCanopy");
    
    active_ = subDict.empty() ? false : true;

    if (active_)
    {
        if (!subDict.found("cellSetName"))
        {
            FatalErrorInFunction << "Must specify a cell set over which to apply wind farm"
                                 << "canopy force, 'cellSetName'."
                                 << abort(FatalError);
        }
        if (!subDict.found("R"))
        {
            FatalErrorInFunction << "Must specify radius from center of hurricane, 'R'."
                                     << abort(FatalError);
        }

       
        cellSetName_ = subDict.lookupOrDefault<word>("cellSetName","windFarmCanopyCellSet");
        cellSet_ = cellSet(mesh_,cellSetName_).toc();

        CtTable_(Function1<scalar>::New("CtTable",subDict));

        R = subDict.lookupOrDefault<scalar>("R",40.0E3);

        Info << "  -using values:" << endl;
        Info << "      R = " << R/1000.0 << " km" << endl;
    }

    else
    {
        Info << "  -no wind-farm canopy specified. Skipping..." << endl;   
    }

}


void Foam::windFarmCanopy::update()
{
    if (active_)
    {

/*
        // Update the source term.
        for (int i = 0; i < nLevels; i++)
        {
            vector sourceAtLevel = Zero;

            scalar V = V_surface + (planeHeights[i] * dVdz);
            scalar dVdR = dVdR_surface + (planeHeights[i] * dVdRdz);

            sourceAtLevel.x() =  Foam::sqr(Ubar[i].x()) / R
                                +Ubar[i].y() * (V/R)
                                -(f*V + (Foam::sqr(V) / R));
            sourceAtLevel.y() = -Ubar[i].x() * dVdR
                                -Ubar[i].x() * (V/R);

          //Info << planeHeights[i] << tab << V_surface << tab << dVdR_surface << tab << V << tab << dVdR << tab << sourceAtLevel << endl;

            forAll(cellsInPlane[i], j)
            {
                source_[cellsInPlane[i][j]] = sourceAtLevel;
            }
        }
*/
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windFarmCanopy::windFarmCanopy
(
    const IOdictionary& dict,
    const volVectorField& U
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    dict_(dict),

    // Set the pointer to the velocity field
    U_(U),

    // Initialize the body force field
    source_
    (
        IOobject
        (
            "windFarmCanopySource",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce",dimVelocity/dimTime,vector::zero)
    ),

    active_(false)
{
    read();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::windFarmCanopy::~windFarmCanopy()
{}


// ************************************************************************* //
