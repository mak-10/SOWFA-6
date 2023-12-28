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

       
        cellSetName_ = subDict.lookupOrDefault<word>("cellSetName","windFarmCanopyCellSet");
        cellSet_ = cellSet(mesh_,cellSetName_).toc();

       // CtTable_(Function1<scalar>::New("CtTable",subDict));

        Cft_ = subDict.lookupOrDefault<scalar>("Cft",0.074);
	canopyHeight_ = subDict.lookupOrDefault<scalar>("canopyHeight",100);

	Info << "  -using values:" << endl;
        Info << "      Cft = " << Cft_ << endl;
	Info << "      canopy Height = " << canopyHeight_ << endl;
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
     	scalar Tc = 0;
	scalar Vc = 0;

	forAll(cellSet_, i)
	{
    		const label cellID = cellSet_[i];

    		scalar T = 0.5 * sqr(mag(U_[cellID])) * Cft_ * (1 / canopyHeight_);
    		source_[cellID] -= T * (U_[cellID] / mag(U_[cellID]));
    		Tc += 1.225 * T * Vcells_[cellID];
		Vc += Vcells_[cellID];
	}

	source_.correctBoundaryConditions();

	reduce(Tc, sumOp<scalar>());
	Info << "Total Thrust in Canopy" << Tc << endl;
	Info << "Total Volume of Canopy" << Vc << endl;

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

    active_(false),
    // Initialize the body force field
    source_
    (
        IOobject
        (
            "windFarmCanopySource",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("bodyForce", dimVelocity/dimTime, vector::zero)
    )

    //Vcells_(
    //    IOobject
    //    (
    //        "Vcells",
    //        runTime_.timeName(),
    //        mesh_,
    //        IOobject::MUST_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    mesh_,
    //    dimensionedScalar("cellVolume", dimVolume, 0.0)
   // )

    //active_(false)
{
    read();
    update();

    // Extract scalar field from mesh_.V() and assign it to Vcells_
    Vcells_ = mesh_.V();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::windFarmCanopy::~windFarmCanopy()
{}


// ************************************************************************* //
