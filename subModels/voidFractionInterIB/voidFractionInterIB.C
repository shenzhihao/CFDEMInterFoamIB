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

#include "voidFractionInterIB.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "dataExchangeModel.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(voidFractionInterIB, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    voidFractionInterIB,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
voidFractionInterIB::voidFractionInterIB
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    alphaLimited_(0),
    scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
    sqrtThree_(sqrt(3.0)),
    vertexSolidFrac(8*sm.mesh().nCells(),0.0)
{
    Info << "\n\n W A R N I N G - do not use in combination with differentialRegion model! \n\n" << endl;
    //Info << "\n\n W A R N I N G - this model does not yet work properly! \n\n" << endl;
    maxCellsPerParticle_=readLabel(propsDict_.lookup("maxCellsPerParticle"));
    //particleCloud_.setMaxCellsPerParticle(readLabel(propsDict_.lookup("maxCellsPerParticle"))); // alternative to line above

    if(scaleUpVol_ < 1){ FatalError<< "scaleUpVol shloud be > 1."<< abort(FatalError); }
    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin shloud be > 1 and < 0.01." << abort(FatalError); }
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voidFractionInterIB::~voidFractionInterIB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void voidFractionInterIB::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes,double**& particleV) const
{

    int numprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    const boundBox& globalBb = particleCloud_.mesh().bounds();

    reAllocArrays();

    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);

    vertexSolidFrac=scalarField(8*particleCloud_.mesh().nCells(),0.0);

    for(int index=0; index < particleCloud_.numberOfParticles(); index++)
    {
        // 遍历particlecloud中的所有颗粒点
        //if(mask[index][0])
        //{
            //reset
            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0;
                particleVolumes[index][subcell] = 0;
            }
            cellsPerParticle_[index][0]=1.0;
            particleV[index][0]=0;

            //collecting data 收集颗粒信息
            label particleCenterCellID = particleCloud_.cellIDs()[index][0];
            // 获取颗粒中心所在的网格单元编号particleCenterCellID。通过locateModel查找获得
            scalar radius = particleCloud_.radius(index);
            // 获得颗粒半径radius
            vector positionCenter = particleCloud_.position(index);
            // 获得颗粒中心三维位置positionCenter

            if (particleCenterCellID >= 0)
            {
                labelHashSet hashSett;

                //compute the voidfraction for the cell "particleCentreCellID
                vector cellCentrePosition = particleCloud_.mesh().C()[particleCenterCellID];
                // 获得颗粒中心所在的网格单元的中心坐标cellCentrePosition

                scalar fcell = pointInParticle(index, positionCenter, cellCentrePosition);
                // 获得颗粒中心和网格中心距离比上颗粒index的半径的比值的平方 再减去一
                vector minPeriodicParticlePos=positionCenter;
                // 颗粒中心坐标
                if(particleCloud_.checkPeriodicCells()) //consider minimal distance to all periodic images of this particle
                {
                    fcell = minPeriodicDistance(index,cellCentrePosition, positionCenter, globalBb,
                                             minPeriodicParticlePos,
                                             particleCloud_.wall_periodicityCheckRange());
                }
                scalar centreDist=mag(cellCentrePosition-minPeriodicParticlePos);
                // 颗粒中心与网格单元中心距离

                scalar corona = 0.5*sqrtThree_*cbrt(particleCloud_.mesh().V()[particleCenterCellID]);
                // 对于正方体，计算网格中心到网格节点的距离，网格对角线的一半。非正方体计算值小于对角线一半
                vector coronaPoint = cellCentrePosition;
                // 网格中心点
                if(centreDist > 0.0)
                {
                  // 颗粒与网格单元中心不重合，大部分情况不会重合
                  coronaPoint = cellCentrePosition + (cellCentrePosition - minPeriodicParticlePos) * (corona / centreDist);
                  // 计算颗粒和网格单元中心连线与网格单元包络球体的交点，并重新赋值给coronaPoint
                }

                if(pointInParticle(index, minPeriodicParticlePos, coronaPoint) < 0.0)
                {
                  // coronaPoint被颗粒包围，说明网格单元被完全包围
                  voidfractionNext_[particleCenterCellID] = 0;
                }
                else
                {
                  // coronaPoint没有被颗粒包围，说明网格单元可能没有被完全包围
                    const labelList& vertices = particleCloud_.mesh().cellPoints()[particleCenterCellID];
                    // 获取单元的顶点序列
                    // int nn = 0.0;
                    // forAll(vertices, i) nn ++;

                    scalar ratio = 0.125; //1.0 / static_cast<double>(nn);
                    forAll(vertices, i)
                    {
                        vector vertexPosition = particleCloud_.mesh().points()[vertices[i]];
                        // 节点坐标
                        if (i>7) i=7;
                        scalar oldSolidFracOnVertex=vertexSolidFrac[8*particleCenterCellID+i];
                        // 获取记录值

                        scalar fvertex = pointInParticle(index, positionCenter, vertexPosition);
                        // 判断网格节点与颗粒相对位置，fvertex<0表示节点被颗粒包围，否这没有被包围
                        if(particleCloud_.checkPeriodicCells()) 
                        { //consider minimal distance to all periodic images of this particle
                            fvertex = minPeriodicDistance(index,vertexPosition, positionCenter, globalBb, 
                                                     minPeriodicParticlePos,
                                                     particleCloud_.wall_periodicityCheckRange());
                        }

                        if(fcell < 0.0 &&  fvertex < 0.0)
                        {
                          // 网格中心与节点中心都被颗粒包围
                          // 贡献值为ratio
                            if(ratio>oldSolidFracOnVertex)
                            {
                              voidfractionNext_[particleCenterCellID]+=oldSolidFracOnVertex;
                              voidfractionNext_[particleCenterCellID]-=ratio;
                              vertexSolidFrac.replace(8*particleCenterCellID+i,ratio);
                            }
                        }
                        else if(fcell < 0.0 && fvertex >= 0.0) 
                        {
                          // 网格中心被颗粒包围，节点中心没有被颗粒包围
                          // 贡献值为ratio*lambda
                            //compute lambda
                            scalar lambda = segmentParticleIntersection(index, minPeriodicParticlePos, cellCentrePosition, vertexPosition);
                            // 计算从网格中心开始的长度
                            if(ratio*lambda > oldSolidFracOnVertex)
                            {
                              voidfractionNext_[particleCenterCellID]+=oldSolidFracOnVertex;
                              voidfractionNext_[particleCenterCellID]-=ratio*lambda;
                              vertexSolidFrac.replace(8*particleCenterCellID+i,ratio*lambda);
                            }
                            // voidfractionNext_[particleCenterCellID] -= ratio*lambda;
                        } else if(fcell >= 0.0 && fvertex < 0.0) 
                        {
                          // 网格中心没有被颗粒包围，节点中心被颗粒包围
                            //compute lambda
                            scalar lambda = segmentParticleIntersection(index, minPeriodicParticlePos, vertexPosition, cellCentrePosition);
                            // 计算从节点开始的长度
                            if(ratio*lambda > oldSolidFracOnVertex)
                            {
                              voidfractionNext_[particleCenterCellID]+=oldSolidFracOnVertex;
                              voidfractionNext_[particleCenterCellID]-=ratio*lambda;
                              vertexSolidFrac.replace(8*particleCenterCellID+i,ratio*lambda);
                            }
                            // voidfractionNext_[particleCenterCellID] -= ratio*lambda;
                        }
                    }
                } //end particle partially overlapping with cell

                //generating list with cell and subcells
                // 递归调用，从中心单元cell出发递归地向四周扩散处理
                buildLabelHashSet(index,
                                  minPeriodicParticlePos, 
                                  particleCenterCellID, 
                                  hashSett, 
                                  vertexSolidFrac,
                                  true);

                //Add cells of periodic particle images on same processor
                if(particleCloud_.checkPeriodicCells())
                {
                  int doPeriodicImage[3];
                  for(int iDir=0;iDir<3;iDir++)
                  {
                    doPeriodicImage[iDir]= 0;
                    if( (minPeriodicParticlePos[iDir]+radius)>globalBb.max()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0)
                    {
                       doPeriodicImage[iDir] =-1;
                    }
                    if( (minPeriodicParticlePos[iDir]-radius)<globalBb.min()[iDir] && particleCloud_.wall_periodicityCheckRange(iDir)>0)
                    {
                         doPeriodicImage[iDir] = 1;
                    }
                  }
                  //scan directions and map particles
                  List<vector> particlePosList;         //List of particle center position
                  List<label>  particleLabelList;

                  int copyCounter=0;
                  // Note: for other than ext one could use xx.append(x)
                  // instead of setSize
                  particlePosList.setSize(particlePosList.size()+1, minPeriodicParticlePos);

                  //x-direction
                  if(doPeriodicImage[0]!=0) {
                    particlePosList.setSize(particlePosList.size()+1, particlePosList[copyCounter]
                                              + vector(
                                                               static_cast<double>(doPeriodicImage[0])
                                                              *(globalBb.max()[0]-globalBb.min()[0]),
                                                              0.0,
                                                              0.0)
                                               );
                    copyCounter++;
                  }
                  //y-direction
                  int currCopyCounter=copyCounter;
                  if(doPeriodicImage[1]!=0) {
                    for(int yDirCop=0; yDirCop<=currCopyCounter; yDirCop++) {
                      particlePosList.setSize(particlePosList.size()+1, particlePosList[yDirCop]
                          + vector(
                                                              0.0,
                                                               static_cast<double>(doPeriodicImage[1])
                                                              *(globalBb.max()[1]-globalBb.min()[1]),
                                                              0.0)
                                               );
                      copyCounter++;
                     }
                  }
                  //z-direction
                  currCopyCounter=copyCounter;
                  if(doPeriodicImage[2]!=0) {
                    for(int zDirCop=0; zDirCop<=currCopyCounter; zDirCop++) {
                      particlePosList.setSize(particlePosList.size()+1, particlePosList[zDirCop]
                          + vector(
                                                              0.0,
                                                              0.0,
                                                               static_cast<double>(doPeriodicImage[2])
                                                              *(globalBb.max()[2]-globalBb.min()[2])
                                                              )
                                               );
                       copyCounter++;
                     }
                  }

                  //add the nearest cell labels
                  particleLabelList.setSize(particleLabelList.size()+1,particleCenterCellID);
                  for(int iPeriodicImage=1;iPeriodicImage<=copyCounter; iPeriodicImage++)
                  {
                    label partCellId =

                        particleCloud_.mesh().findNearestCell(particlePosList[iPeriodicImage]);
                    particleLabelList.setSize(particleLabelList.size()+1,partCellId);

                    buildLabelHashSet(index,
                                      particlePosList[iPeriodicImage], 
                                      particleLabelList[iPeriodicImage], 
                                      hashSett, 
                                      vertexSolidFrac,
                                      false);
                  }
                } //end particleCloud_.checkPeriodicCells()

                scalar hashSetLength = hashSett.size();
                if (hashSetLength > maxCellsPerParticle_) {
                    FatalError<< "big particle algo found more cells ("<< hashSetLength 
                              <<") than storage is prepared ("<<maxCellsPerParticle_<<")" << abort(FatalError);
                } else if (hashSetLength > 0) {
                  cellsPerParticle_[index][0]=hashSetLength;
                  hashSett.erase(particleCenterCellID);

                  for(label i=0;i<hashSetLength-1;i++) {
                    label cellI=hashSett.toc()[i];
                    particleCloud_.cellIDs()[index][i+1]=cellI; //adding subcell represenation
                  }
                }//end cells found on this proc
            }// end found cells
        //}// end if masked
    }// end loop all particles

    for(label index=0; index < particleCloud_.numberOfParticles(); index++)
    {
      for(label subcell=0;subcell<cellsPerParticle_[index][0];subcell++) 
      {
        label cellID = particleCloud_.cellIDs()[index][subcell];

        if(cellID >= 0) {
          if(voidfractionNext_[cellID] < 0.0)
            voidfractionNext_[cellID] = 0.0;
          voidfractions[index][subcell] = voidfractionNext_[cellID];
        } else {
          voidfractions[index][subcell] = -1.;
        }
      }
    }
}

void voidFractionInterIB::buildLabelHashSet
(
    int index,
    const vector position,
    const label cellID,
    labelHashSet& hashSett, 
    scalarField& vertexsolidfrac,
    bool initialInsert //initial insertion of own cell
)const
{
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    if(initialInsert) hashSett.insert(cellID);
    // 首先将当前网格ID加进到hashSet之中去

    const labelList& neighborcells = particleCloud_.mesh().cellCells()[cellID];

    forAll(neighborcells,i) 
    {
        scalar Solidfrac = 0.0;

        label neighbor=neighborcells[i];
        vector cellCentrePosition = particleCloud_.mesh().C()[neighbor];
        scalar centreDist = mag(cellCentrePosition-position);
        scalar fcell = pointInParticle(index, position, cellCentrePosition);
        scalar corona = 0.5*sqrtThree_*cbrt(particleCloud_.mesh().V()[neighbor]);

        vector coronaPoint = cellCentrePosition;
        if(centreDist > 0.0)
          coronaPoint = cellCentrePosition + (cellCentrePosition - position) * (corona / centreDist);
        

        if(!hashSett.found(neighbor) && pointInParticle(index, position, coronaPoint) < 0.0)
        {
          // !hashSett.found(neighbor)表示neighbor没有被加进到hashset中过，即尚未被处理
          // pointInParticle<0表示coronoint在颗粒内，网格单元被颗粒完全包围，voidfraction直接赋值为0
            voidfractionNext_[neighbor] = 0;
            
            buildLabelHashSet(index, 
                              position, 
                              neighbor, 
                              hashSett, 
                              vertexsolidfrac,
                              true);
        }
        else if(!hashSett.found(neighbor)) 
        {
          // !hashSett.found(neighbor)表示neighbor没有被加进到hashset中过，即尚未被处理
          // pointInParticle>=0表示coronoint不在颗粒内，网格单元没有被颗粒完全包围，voidfraction需要计算
            // int nn = 0.0;
            // forAll(vertexPoints, i) {
            //   nn ++;
            // }

            const labelList& vertexPoints = particleCloud_.mesh().cellPoints()[neighbor];
            scalar ratio = 0.125; //1.0 / static_cast<double>(nn);
            forAll(vertexPoints, j) 
            {
                if(j>7) j=7;

                vector vertexPosition = particleCloud_.mesh().points()[vertexPoints[j]];

                // 获取solidFrac记录值

                scalar fvertex = pointInParticle(index, position, vertexPosition);

                if (fcell < 0.0 && fvertex < 0.0) 
                {
                  // 网格中心与节点中心都被颗粒包围
                  // 贡献值为ratio
                  Solidfrac += ratio;
                }
                else if (fcell < 0.0 && fvertex > 0.0)
                {
                  //网格单元中心被包围，节点没有被包围
                  scalar lambda = segmentParticleIntersection(index, position, cellCentrePosition, vertexPosition);
                  Solidfrac += ratio * lambda;
                }
                else if (fcell > 0.0 && fvertex < 0.0) 
                {
                  //网格单元中心没有被包围，节点被包围
                  scalar lambda = segmentParticleIntersection(index, position, vertexPosition, cellCentrePosition);
                  Solidfrac += ratio * lambda;
                }
            }

            if(Solidfrac > 1.0)
              Solidfrac = 1.0;

            scalar newViod = 1.0 - Solidfrac;
            if(newViod < voidfractionNext_[neighbor])
            voidfractionNext_[neighbor] = newViod;

            // 限制下限为0
            if(voidfractionNext_[neighbor] < 0) 
              voidfractionNext_[neighbor] = 0;

            // if(!(Solidfrac == 0.0))
            if(!(voidfractionNext_[neighbor] == 1.0))
            {
              // !(Solidfrac == 0.0)表示网格的neighbor可能与颗粒相交需要进一步处理
              buildLabelHashSet(index, 
                                position, 
                                neighbor, 
                                hashSett, 
                                vertexsolidfrac,
                                true);
            }
        }
    }
}

double voidFractionInterIB::segmentParticleIntersection(int index, vector positionCenter, vector pointInside, vector pointOutside) const
{
  scalar radius =  particleCloud_.radius(index);
  scalar a = (pointOutside - pointInside)&(pointOutside - pointInside);
  scalar b = 2.*(pointOutside - pointInside)&(pointInside - positionCenter);
  scalar c = ((pointInside - positionCenter)&(pointInside - positionCenter)) - radius*radius;
  scalar lambda_ = 0.0;
  scalar lambda = 0.0;
  scalar D = b*b - 4.0*a*c;
  double eps = 1e-12;
  if(D >= 0.0) {
    double sqrtD = sqrt(D);
    lambda_ = (-b + sqrtD)/(2.0*a);
    if(lambda_ >= -eps && lambda_ <= 1.0+eps)
      lambda = lambda_;
    else {
      lambda_ = (-b - sqrtD)/(2.0*a);
      if(lambda_ >= -eps && lambda_ <= 1.0+eps)
        lambda = lambda_;
    }
  }

  if(lambda < 0.0)
    lambda = 0.0;
  if(lambda > 1.0)
    lambda = 1.0;
  return lambda;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
