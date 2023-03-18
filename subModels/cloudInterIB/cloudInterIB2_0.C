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

#include "fileName.H"
#include "cloudInterIB2_0.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "mpi.h"
#include "IOmanip.H"
#include "OFversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cloudInterIB20::cloudInterIB20
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    angularVelocities_(NULL),
    DEMTorques_(NULL),
    pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
    pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))),
    haveEvolvedOnce_(false),
    skipLagrangeToEulerMapping_(false),
    skipAfter_(false),
    timeStepsToSkip_(0),
    calculateTortuosity_(false),
    frontMeshRefine_(false),
    IBDragPerV_
    (
        IOobject
        (
            "IBDragPerV_cloud",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("IBDragPerV_cloud", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    IBDragPresPerV_
    (
        IOobject
        (
            "IBDragPresPerV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("IBDragPresPerV", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    ),
    IBDragViscPerV_
    (
        IOobject
        (
            "IBDragViscPerV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("IBDragViscPerV", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    )
{
    Info << "\n==========Start Constructor of myCloudinterIB2.0============\n" << endl;
    if(this->couplingProperties().found("skipLagrangeToEulerMapping"))
    {
        Info << "Will skip lagrange-to-Euler mapping..." << endl;
        skipLagrangeToEulerMapping_=true;
    }
    if(this->couplingProperties().found("timeStepsBeforeSkipping"))
    {
        skipAfter_=true;
        timeStepsToSkip_ =  readScalar
        (
            this->couplingProperties().lookup("timeStepsBeforeSkipping")
        );
        Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }
    if(this->couplingProperties().found("tortuosity"))
    {
        calculateTortuosity_ = true;
        flowDir_ = this->couplingProperties().subDict("tortuosity").lookup("flowDirection");
        flowDir_ = flowDir_ / mag(flowDir_);
        Info << "Will calculate tortuosity in the mean flow direction ("<<flowDir_[0]<<" "<<flowDir_[1]<<" "<<flowDir_[2]<<")"<< endl;
    }

    //Must check for walls in case of checkPeriodicCells
    //periodic check will mirror particles and probing points to ensure proper behavior near processor bounds
    if(checkPeriodicCells_)
    {
    //Enforce reading of the blocking for periodic checks
    if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("x")))
        wall_periodicityCheckRange_[0] = 0;
    if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("y")))
        wall_periodicityCheckRange_[1] = 0;
    if(readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("z")))
        wall_periodicityCheckRange_[2] = 0;
    
    if(this->couplingProperties().found("wall_periodicityCheckTolerance"))
        wall_periodicityCheckTolerance_ = readScalar (this->couplingProperties().lookup("wall_periodicityCheckTolerance"));
    }
    Info << "\nEnd Constructor of myCloudinterIB2.0\n" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cloudInterIB20::~cloudInterIB20()
{
    dataExchangeM().destroy(angularVelocities_,3);
    dataExchangeM().destroy(dragPrev_,3);
    dataExchangeM().destroy(DEMTorques_,3);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::cloudInterIB20::getDEMdata()
{
    cfdemCloud::getDEMdata();
    dataExchangeM().getData("omega","vector-atom",angularVelocities_);
}

bool Foam::cloudInterIB20::reAllocArrays() const
{
    if(cfdemCloud::reAllocArrays())
    {
        Info <<"Foam::cloudInterIB20::reAllocArrays()"<<endl;
        dataExchangeM().allocateArray(angularVelocities_,0,3);
        dataExchangeM().allocateArray(dragPrev_,0,3);
        dataExchangeM().allocateArray(DEMTorques_,0,3);
        return true;
    }
    return false;
}

void Foam::cloudInterIB20::giveDEMdata()
{

    cfdemCloud::giveDEMdata();
    dataExchangeM().giveData("hdtorque","vector-atom", DEMTorques_);
}

inline double ** Foam::cloudInterIB20::DEMTorques() const
{
    return DEMTorques_;
}


bool Foam::cloudInterIB20::evolve
(
    volScalarField& voidfraction,
    volScalarField& alpha1,
    volScalarField& interFace,
    volVectorField& ibdragperv,
    volVectorField& ibdragpressperv,
    volVectorField& ibdragviscoperv,
    volVectorField& U,
    volScalarField& p
)
{
    Info << "===============================================" << endl;
    Info << "Start myCloudinterib2.0 evolve" << endl;
    numberOfParticlesChanged_ = false;
    arraysReallocated_=false;
    bool doCouple=false;

    // Info << "skipAfter = " << skipAfter_ << endl;
    if(skipAfter_) {
      if(timeStepsToSkip_<1)
        skipLagrangeToEulerMapping_=true;
    }

    if(!writeTimePassed_ && mesh_.time().outputTime())
        writeTimePassed_=true;

    // doCoupleNow()计算当前时间步应不应该进行耦合
    if (dataExchangeM().doCoupleNow())
    {
        dataExchangeM().couple(0);
        doCouple=true;
        if(!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_)
        {
            // 获取颗粒坐标，半径，速度，id等信息，具体取决于dataexchagemodel
            getDEMdata();
            // 对颗粒中心所在单元进行定位，获取其单元ID。定位算法是OF的meshsearch.H,八叉树octree
            locateM().findCell(NULL,positions_,cellIDs_,numberOfParticles());
            Info <<"locate model done\n";
            // 计算voidfraction
            voidFractionM().setvoidFraction(NULL,voidfractions_,particleWeights_,particleVolumes_,particleV_);
            // printf("set void fraction done\n");
            // 设置用于动网格更新的interface变量
            setInterFace(interFace,alpha1);
        }
        // update voidFractionField更新voidfraction场
        // voidfractions_是voidFractionModel中的值，而voidfraction是此函数的实参，表示主程序中的voidfraction
        voidfraction == voidFractionM().voidFractionNext(); // there might be a better approach, see cfdemCloud.C
        voidfraction.correctBoundaryConditions();

        // set particles forces
        // 首先对变量初始化为0，清除上一步的计算结果
        for(int index = 0;index <  numberOfParticles_; ++index)
        {
            for(int i=0;i<3;i++)
            {
                impForces_[index][i] = 0;
                expForces_[index][i] = 0;
                DEMForces_[index][i] = 0;
            }
        }
        for (int i=0;i<nrForceModels();i++)
        {
            // 通过forcemodel设置颗粒所受的力
            forceM(i).setForce();
        }
        
        // 获得ibdrag信息
        // 此处用于数据输出，可注释掉
        Info << "Get ibdragperv" << endl;
        ibdragperv = interIBDragPerV(U,p);
        ibdragpressperv = IBDragPressPerV(U,p);
        ibdragviscoperv = IBDragViscoPerV(U,p);
        // 将计算得到的力传回DEM程序部分
        giveDEMdata();

        dataExchangeM().couple(1);      
        haveEvolvedOnce_=true;
    }
    // do particle IO
    Info << "Dump DEM data" << endl;
    IOM().dumpDEMdata();

    if(skipAfter_)
    {
        timeStepsToSkip_--;
        Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }
    Info << "Finish myCloudinterib2.0 evolve" << endl;
    Info << "__________________________________________________" << endl;
    return doCouple;
}

//defines the mesh refinement zone around a particle
//twice the particle size in each direction
void Foam::cloudInterIB20::setInterFace
(
    volScalarField& interFace,
    volScalarField& alpha1
)
{
    interFace == dimensionedScalar("zero", interFace.dimensions(), 0.);
    for(int par=0; par< numberOfParticles(); par++)
    {
        vector ParPos(positions()[par][0],positions()[par][1],positions()[par][2]);
        const boundBox& globalBb = mesh().bounds();
        double skin = 2.0;
        forAll(mesh_.C(),cellI)
        {
            vector posC = mesh_.C()[cellI];
            if(checkPeriodicCells_)
            {
                // Some cells may be located on the other side of a periodic boundary.
                // In this case, the particle center has to be mirrored in order to correctly
                // evaluate the interpolation points.
                vector minPeriodicParticlePos=ParPos;
                voidFractionM().minPeriodicDistance(par,posC, ParPos, globalBb, 
                                                    minPeriodicParticlePos, 
                                                    wall_periodicityCheckRange());

                ParPos = minPeriodicParticlePos;
            }
            double value = voidFractionM().pointInParticle(par, ParPos, posC, skin);
            if(value <= 0.0)
            {
                interFace[cellI] = value + 1.0;
            }
            if(alpha1[cellI]  > 0.2 && alpha1[cellI] < 0.8)
            {
                interFace[cellI] = max(alpha1[cellI],interFace[cellI]);
            }
        }
    }
}

void Foam::cloudInterIB20::calcVelocityCorrection
(
    volScalarField& rho,
    volScalarField& p,
    volVectorField& U,
    volScalarField& phiIB,
    volScalarField& voidfraction,
    volScalarField& udivmid,
    volVectorField& usdoamin,
    volScalarField& cella,
    volScalarField& alpha1
)
{
    Info << "__________________________________________________" << endl;
    Info << "start velocity correction" << endl;
    setParticleVelocity(
        U,
        usdoamin,
        cella,
        alpha1);

    udivmid = fvc::div(U);
    scalar doDivCor = couplingProperties_.lookupOrDefault<scalar>("doDivCor",0.);
    // scalar 必须用小数0.,不能用整数0
    // dictionary.lookup("a")|lookupOrDefault("a","b")
       
    if(doDivCor == 1)
    {
        Info << "doConvCorr == 1, phiIBEqn" << endl;
        fvScalarMatrix phiIBEqn
        (
            fvm::laplacian(phiIB) == fvc::div(U)
        );
        if(phiIB.needReference()) 
        {
             phiIBEqn.setReference(pRefCell_, pRefValue_);
        }  
        phiIBEqn.solve();    
        // U = fvc::reconstruct
        // (fvc::flux(U) - fvc::snGrad(phiIB)*mesh_.magSf());
        U = U - fvc::grad(phiIB);        
        U.correctBoundaryConditions();
        p = p + rho * phiIB/U.mesh().time().deltaT();
        p.correctBoundaryConditions();

        Info << "End correct p and U" << endl;
    }    
//---------------------------------------------------------------
//checkinterface
    // if (couplingProperties_.found("checkinterface"))
    // {
    //     Info << "checking no-slip on interface..." << endl;
    // }
    // Info << "stop velocity correction" << endl;
//---------------------------------------------------------------
    Info << "End velocity correction" << endl;
    Info << "__________________________________________________" << endl;
}

// 将DEM计算结果映射到CFD网格中
// 设置权重 ???
void Foam::cloudInterIB20::setParticleVelocity
(
    volVectorField& U,
    volVectorField& usdomain,
    volScalarField& cella,
    volScalarField& alpha1
)
{
    Info << "__________________________________________________" << endl;
    Info << "setParticleVelocity" << endl;
    label cellI = 0;
    bool useCoeV = false;
    vector uParticle(0,0,0);
    vector rVec(0,0,0);
    vector velRot(0,0,0);
    vector angVel(0,0,0);

    // double correff = 0.5;
    scalar CoeV_local = readScalar(couplingProperties_.lookup("Coe_V_local"));
    scalar CoeV_global = readScalar(couplingProperties_.lookup("Coe_V_global"));
    scalar real_CoeV;

    cella = 0;
    usdomain = dimensionedVector("usdomain", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0));

    Info << "Total number of particles = " << numberOfParticles() << endl;

    for(int index=0; index < numberOfParticles(); index++)
    // index最大值表示总颗粒数
    {
        useCoeV = false;
        for(int subCell=0; subCell < cellsPerParticle()[index][0]; subCell++)
        {
            cellI = cellIDs()[index][subCell];
            if (cellI >= 0)
            {
                if (alpha1[cellI] > 0.02 && alpha1[cellI] < 0.98)
                {
                    useCoeV = true;
                    break;
                }
            }
        }

        if (useCoeV)
            real_CoeV = CoeV_local;
        else
            real_CoeV = CoeV_global; 

        // Info << "\nTotal subcells for particle " << index << " is " << cellsPerParticle()[index][0] << endl;
        for(int subCell=0; subCell < cellsPerParticle()[index][0]; subCell++)
        // subCell最大值表示每个颗粒所占据的CFD网格数
        {
            //Info << "subCell=" << subCell << endl;
            cellI = cellIDs()[index][subCell];
            // Info << "\nThe cell ID of subCell " << subCell << 
            // " is " << cellI << endl;
            if (cellI >= 0)
            {
                if (!cella[cellI])
                {
                    cella[cellI] = 1;
                    // calc particle velocity
                    for(int i=0;i<3;i++) rVec[i] = U.mesh().C()[cellI][i] - position(index)[i];
                    for(int i=0;i<3;i++) angVel[i] = angularVelocities()[index][i];
                    velRot = angVel^rVec;
                    for(int i=0;i<3;i++) uParticle[i] = velocities()[index][i] + velRot[i];
    
                    // impose field velocity
                    // (1-voidfractions_[index][subCell])*uParticle+
                    // voidfractions_[index][subCell]*U[cellI];
                    
                    usdomain[cellI] = 
                    real_CoeV * 
                    (1-voidfractions_[index][subCell]) *
                    U[cellI]
                    ;
    
                    U[cellI] = 
                    real_CoeV * 
                    (1-voidfractions_[index][subCell]) *
                    uParticle +
                    voidfractions_[index][subCell] *
                    U[cellI]
                    ;
                }            
            }
        }
    }
    U.correctBoundaryConditions();
    Info << "End of setParticleVelocity" << endl;  
    Info << "__________________________________________________" << endl;
}

vector Foam::cloudInterIB20::angularVelocity(int index)
{
    vector vel;
    for(int i=0;i<3;i++) vel[i] = angularVelocities_[index][i]; 
    return vel;
}

double Foam::cloudInterIB20::getTortuosity(vector dir)
{
    volVectorField U = mesh_.lookupObject<volVectorField>("U");
    volScalarField voidfraction = mesh_.lookupObject<volScalarField>("voidfraction");
    double ux = 0.0;
    double umag = 0.0;
    forAll(mesh_.V(),cellI)
    {
        if(voidfraction[cellI] > 0.5)
        {
            double V = mesh_.V()[cellI];
            ux += ((U[cellI] & dir))*V;
            umag += mag(U[cellI])*V;
        }
    }
    //double ux_reduced = 0.0;
    //double umag_reduced = 0.0;
    //MPI_Allreduce(&ux, &ux_reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&umag, &umag_reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    reduce(umag, sumOp<scalar>());
    reduce(ux, sumOp<scalar>());
    double tortuosity = ux == 0.0 ? 1.0 : umag / ux;
    return tortuosity;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void Foam::cloudInterIB20::setRefinementField(volScalarField* refine_)
{
 //Function to allow for setting and activating special refinement operations
 frontMeshRefineField_ = refine_;
 frontMeshRefine_ = true;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

///////////////////////////////////////////////////////////////////////////////////////////////////////

const volVectorField& Foam::cloudInterIB20::interIBDragPerV(const volVectorField& U,const volScalarField& p)
{
    Info << "In revised interIBdragPerV" << endl;
    #ifdef compre
        IBDragPerV_ = forceM(0).forceSubM(0).muField()*fvc::laplacian(U)-fvc::grad(p);
    #else
        // Info << "laplacian(mu(),U)" << endl;
        IBDragPerV_ = fvc::laplacian(forceM(0).forceSubM(0).rhoField()*forceM(0).forceSubM(0).nuField(),U);
        // Info << "Grad(p)" << endl;
        IBDragPerV_ += -fvc::grad(p);
    #endif
    Info << "Out revised interIBdragPerV" << endl;
    return IBDragPerV_;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

const volVectorField& Foam::cloudInterIB20::IBDragPressPerV(const volVectorField& U,const volScalarField& p)
{
    // Info << "Get ibdrag pressure part" << endl;
    Info << "IBDragPressPerV = grad(p)" << endl;
    IBDragPresPerV_ = -fvc::grad(p);
    // Info << "Out IBDragPressPerV" << endl;
    return IBDragPresPerV_;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////

const volVectorField& Foam::cloudInterIB20::IBDragViscoPerV(const volVectorField& U,const volScalarField& p)
{
    // Info << "Get ibdrag viscosity part" << endl;

    Info << "IBDragViscoPerV = laplacian(mu(),U)" << endl;
    IBDragViscPerV_ = fvc::laplacian(forceM(0).forceSubM(0).rhoField()*forceM(0).forceSubM(0).nuField(),U);
    
    // Info << "Out IBDragViscoPerV" << endl;
    return IBDragViscPerV_;
}


} // End namespace Foam

// ************************************************************************* //
