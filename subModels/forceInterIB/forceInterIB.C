/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "forceInterIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(forceInterIB, 0);

addToRunTimeSelectionTable
(
    forceModel,
    forceInterIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceInterIB::forceInterIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    twoDimensional_(false),
    depth_(1),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    pressureFieldName_(propsDict_.lookup("pressureFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_)),
    voidFieldName_(propsDict_.lookup("voidFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidFieldName_)),
    // phaseName_(propsDict_.lookup("phaseName")),
    // phase_(sm.mesh().lookupObject<volScalarField> (phaseName_)),
    interIBDragPerV_
    (
        IOobject
        (
            "interIBDragPerV_",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedVector("interIBDragPerV_", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    useTorque_(false)
{
    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, typeName+".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().writeHeader();


    if (propsDict_.found("twoDimensional"))
    {
        twoDimensional_=true;
        depth_ = readScalar(propsDict_.lookup("depth"));
        Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
        Info << "depth of domain is assumed to be :" << depth_ << endl;
    }

    if(propsDict_.found("useTorque")) useTorque_ = true;

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceInterIB::~forceInterIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forceInterIB::setForce() const
{
    Info << "__________________________________________________" << endl;
    Info << "setforce version: forceInterIB 2.0" << endl;
    label cellI;
    vector drag;
    vector torque;

    scalar Exdrag = readScalar(dict_.lookup("Exdrag"));
    scalar dragcorrcoe = readScalar(dict_.lookup("dragcorrcoe"));

    volVectorField h = calcInterIBDragPerV(U_,p_);
    // *voidfraction_
    ;

    #include "setupProbeModel.H"

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            // drag=vector::zero;
            drag    =   vector(0.,0.,0.);
            torque  =   vector(0.,0.,0.);
            vector  positionCenter = particleCloud_.position(index);
            
            if(forceSubM(0).verbose())
                Info << "\nTotalCell for particle " << index
                << " at " << positionCenter
                << " is " << particleCloud_.cellsPerParticle()[index][0] << endl;

            for(int subCell=0;subCell<particleCloud_.cellsPerParticle()[index][0];subCell++)
            {
                cellI = particleCloud_.cellIDs()[index][subCell];
                if (cellI > -1) // particle Found
                {
                	if(forceSubM(0).verbose())
                        Info << "cellid =" << cellI << "||";
                    vector rc = particleCloud_.mesh().C()[cellI];
                    //drag   += h[cellI]*h.mesh().V()[cellI];
                    drag   += h[cellI]*h.mesh().V()[cellI]
                    * (1 - voidfraction_[cellI])
                    * dragcorrcoe
                    ;
                    torque += (rc - positionCenter)^h[cellI]*h.mesh().V()[cellI]
                    * (1 - voidfraction_[cellI])
                    * dragcorrcoe
                    ;
                    // 针对multisphere模型,防止重复计算
                    h[cellI] = vector(0.,0.,0.);

                }
            }
            if(Exdrag == 1)
                Info << "\ndrag on particle " << index 
                << " is " << drag << endl;

            // set force on particle
            if(twoDimensional_) drag /= depth_;

            //Set value fields and write the probe
            if(probeIt_)
            {
                #include "setupProbeModelfields.H"
                // Note: for other than ext one could use vValues.append(x)
                // instead of setSize
                vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                particleCloud_.probeM().writeProbe(index, sValues, vValues);
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,vector::zero);
            // 此函数内部通过drag变量来赋值
                // impForces_
                // expForces_
                // DEMForces_
            // 这三个变量
            // Info << "drag =" << drag << endl;

            // if(forceSubM(0).verbose()) 
            if(Exdrag == 1) 
                Info << "impForces = " 
                <<impForces()[index][0]<<","
                <<impForces()[index][1]<<","
                <<impForces()[index][2] << endl;
//            if(useTorque_) 
            for(int j=0;j<3;j++) 
                particleCloud_.DEMTorques()[index][j] = torque[j];
            // 总是使用弯矩
        //}
    }
    Info << "__________________________________________________" << endl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const volVectorField& forceInterIB::interIBDragPerV() const
{
    return interIBDragPerV_;
}


const volVectorField& forceInterIB::calcInterIBDragPerV(const volVectorField& U,const volScalarField& p) const
{
    #ifdef compre
        interIBDragPerV_ = forceSubM(0).muField()*fvc::laplacian(U)-fvc::grad(p);
    #else
        // Info << "laplacian(mu(),U)" << endl;
        interIBDragPerV_ = fvc::laplacian(forceSubM(0).rhoField()*forceSubM(0).nuField(),U);
        // Info << "Grad(p)" << endl;
        interIBDragPerV_ += -fvc::grad(p);
    #endif
    return interIBDragPerV_;
}

} // End namespace Foam

// ************************************************************************* //
