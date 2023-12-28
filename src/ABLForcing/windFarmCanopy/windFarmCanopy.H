/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::windFarmCanopy

Description
    A simple wind-farm drag canopy-like model. Toposet it used to create a cell
    set of cells over which the drag force will be applied.  A wind-farm thrust
    coefficient is supplied and used to compute a drag force locally.

SourceFiles
    windFarmCanopy.C

\*---------------------------------------------------------------------------*/

#ifndef windFarmCanopy_H
#define windFarmCanopy_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "cellSet.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class windFarmCanopy Declaration
\*---------------------------------------------------------------------------*/

class windFarmCanopy
{
    // Private data

        //- Constants
            //- Runtime pointer
            const Time& runTime_;
    
            //- Mesh pointer
            const fvMesh& mesh_;

            //- Velocity field pointer
            const volVectorField& U_;

        //- Dictionary
        const IOdictionary& dict_;

        //- Body force field of terms specific to wind-farm canopy.
        volVectorField source_;

	// scalar to access the cell volumes of the cellSet or canopy
	scalarField  Vcells_;

        bool active_;

        //- Parameters for the wind-farm canopy.
        word cellSetName_;
        List<label> cellSet_;
        //autoPtr<Function1<scalar> > CtTable_;
        scalar Cft_;
	scalar canopyHeight_;
	scalar rho_;



    // Private Member Functions

        //- Read a subdict
        void read();


public:

    //- Declare name of the class and its debug switch
    ClassName("windFarmCanopy");


    //- Constructors
    windFarmCanopy
    (
        const IOdictionary& dict,
        const volVectorField& U
    );


    //- Destructor
    virtual ~windFarmCanopy();


    //- Update the source term.
    void update();

    //- Return force
    volVectorField& source()
    {
        return source_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

