/*---------------------------------------------------------------------------*\
Description
    Transient solver for incompressible flow.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6, 
    where additional functionality for CFD-DEM coupling using immersed body
    (fictitious domain) method is added.
\*---------------------------------------------------------------------------*/
// from inter
#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
// #include "multiphaseMixture.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "OFversion.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "pimpleControl.H"
#include "cloudInterIB2_0.H"
#include "implicitCouple.H"
#include "averagingModel.H"
#include "voidFractionModel.H"
#include "dynamicFvMesh.H"
#include "cellSet.H"

int main(int argc, char *argv[])
{
    Info << "Solver name: solverInterIB2.0"<< endl;
    // =========================================
    // Not local time stepping (LTS)
    // 不使用LTS
    // 基本设置
    #include "postProcess.H"

    // 设置并打印基础信息
    #include "setRootCase.H"

    // createTime: 注册时间，注册结构的最根部数据，生成Time类型的runTime对象
    #include "createTime.H"

    // 注册网格，依附于时间
    #include "createDynamicFvMesh.H"

    // 注册piso对象，控制计算过程
    // pisoControl control(mesh);
    pimpleControl control(mesh);

    // 创建时间控制
    #include "createTimeControls.H"

    #include "initContinuityErrs.H"

    // 创建运算所以需要的场，统一在createFields中定义
    #include "createFields.H"

    // 创建耦合cloud*****************
    // Info << "==========Create particlecloud" << endl;
    cloudInterIB20 particleCloud(mesh);

    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time Step = " << runTime.timeName() << nl << endl;
//=================================
//动网格生成
//=================================
        // interFace = mag(mesh.lookupObject<volScalarField>("voidfractionNext"));
        particleCloud.setMeshHasUpdatedFlag(mesh.update()); //dyM
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
        // mixture.correct();
//=================================
//DEM 颗粒受力计算,运动更新
// evolve 内容
// 更新 voidfraction
// 设置颗粒受力
//=================================
        Info << "Update particle motion\n" << endl;
        particleCloud.evolve
        (
            voidfraction, 
            alpha1, 
            interFace, 
            ibdrag, 
            ibdragP, 
            ibdragV, 
            U, 
            p
        );
        Info << "Finish update particle motion\n" << endl;
//=================================
//CFD 流体域计算: PISO 算法求解动量方程和连续方程
//=================================
        if(particleCloud.solveFlow())
        {
            // 类似于pimple.loop()
            Info<< "Solve flow at time: " << runTime.timeName() 
            << nl << endl;
//==============================================================
// 求解相方程
//==============================================================
            // #include "alphaControls.H"
            const dictionary& alphaControls = mesh.solverDict(alpha1.name());
            label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
            label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
            // 结束 include "alphaControls.H"
            //==============================================================
            // 相方程更新
            // #include "alphaEqnSubCycle.H"
            // 此处用到的alphaEqnSubCycle位于/applications/solvers/multiphase/VoF
            if (nAlphaSubCycles > 1)
            {
                // 获取时间不长dt
                dimensionedScalar totalDeltaT = runTime.deltaT();
                // 创建rho phi sum场，初始化为"0"，单位使用与rhophi的单位一致
                surfaceScalarField rhoPhiSum
                (
                    IOobject
                    (
                        "rhoPhiSum",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar("0", rhoPhi.dimensions(), 0.0)
                );
                // tmp<volScalarField> trSubDeltaT;
                // if (LTS)
                // {
                //     trSubDeltaT =
                //         fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
                // }
                for
                (
                    subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                    !(++alphaSubCycle).end();
                )
                {
                    #include "alphaEqn.H"
                    rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
                }

                rhoPhi = rhoPhiSum;
            }
            else
            {
                #include "alphaEqn.H"
            }
            rho == alpha1*rho1 + alpha2*rho2;
            mixture.correct();
// 结束 alphaEqnSubCycle.H
//==============================================================
//==============================================================
// 动量预测步 Momentum predictor
//==============================================================
            Info << "Construct UEqn\n" << endl;
            MRF.correctBoundaryVelocity(U);
            //=================================
            // 构建UEqn Matrix，变量为U
            // 对应偏微分方程: d(rho*U)/dt + div(rho*U^U)
            // 第一项瞬态项 fvm::ddt(rho, U): 欧拉全隐式格式;
            // 第二项对流项 fvm::div(rhoPhi, U): 半隐式格式rhophi为上一时间步
            // 第三项MRF
            // 第四项扩散项 turbulence->divDevRhoReff(rho, U): 粘度项，由turbulence处理
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
                + MRF.DDt(rho, U)
                + turbulence->divDevRhoReff(rho, U)
                ==
                fvOptions(rho, U)
            );

            UEqn.relax();

            fvOptions.constrain(UEqn);

            // 判断是否进行动量预测
            if (control.momentumPredictor())
            {
                Info << "Solve UEqn == Source\n" << endl;
                // reconstruct，由面值构建得到体心值
                //reconstruct: create vol value from sur value
                solve
                (
                    UEqn 
                    ==
                    fvc::reconstruct
                    (
                        (
                            mixture.surfaceTensionForce()
                            - ghf*fvc::snGrad(rho)
                            - fvc::snGrad(p_rgh)
                        ) * mesh.magSf()
                    )
                );

                fvOptions.correct(U);
            }
            // END UEqn
//=================================
// pEqn 求解泊松方程步 PISO算法
//=================================
            Info << "PIMPLE loop\n" << endl;
            while (control.correct())
            {
                // ==============================================
                // 构建HbyA以及其面通量phiHbyA
                // rAU表示reverse of A(U),即 1/A
                // volScalarField rAU = 1.0/UEqn.A();
                volScalarField rAU("rAU", 1.0/UEqn.A());
                // 将1/A插值到网格面上得到(1/A)f
                surfaceScalarField rAUf = fvc::interpolate(rAU);
                // HbyA = H/A
                // UEqn离散后变为 A * U = H,H中包含的U为上一时间步的结果,HbyA得到新的时间步结果
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

                // 得到HbyA的面通量phiHbyA
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                    // 动网格流率修正
                    + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi) 
                );
                MRF.makeRelative(phiHbyA);
                adjustPhi(phiHbyA, U, p_rgh);

                // 速度方程中未包含的部分，更新到phiHbyA中
                surfaceScalarField phiSourcebyA
                (
                    
                    (
                        mixture.surfaceTensionForce()
                        - ghf*fvc::snGrad(rho)
                    ) * rAUf * mesh.magSf()
                );
                phiHbyA += phiSourcebyA;

                constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

                // Non-orthogonal pressure corrector loop
                // 非正交修正步骤，若使用正交网格，此步骤仅执行一次
                while (control.correctNonOrthogonal())
                {   
                    // 构建压力泊松方程并求解
                    fvScalarMatrix p_rghEqn
                    (
                        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
                        // + fvc::ddt(voidfraction)
                    );
                    p_rghEqn.setReference(pRefCell, pRefValue);
                    p_rghEqn.solve(mesh.solver(p_rgh.select(control.finalInnerIter())));
                    //p_rgh -> p_rghEqn

                    if (control.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - p_rghEqn.flux();
                        p_rgh.relax();
                        // U = HbyA + rAU*(fvc::reconstruct(phiSourcebyA/rAUf)-fvc::grad(p_rgh));
                        U = HbyA + rAU * fvc::reconstruct((phiSourcebyA - p_rghEqn.flux())/rAUf);
                        // 最后一步，U=HbyA+1/A*()
                        U.correctBoundaryConditions();
                        fvOptions.correct(U);
                    }
                    // =====================================================
                }
                // 输出连续性误差信息
                #include "continuityErrs.H"

                p = p_rgh + rho*gh;
                
                // 可选的 设置参考点
                if (p_rgh.needReference())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pRefValue - getRefCellValue(p, pRefCell)
                    );
                    p_rgh = p - rho*gh;
                }  
            }
            // 结束PISO算法的求解压力p的循环
            // ddtvoid = fvc::ddt(voidfraction);
        } 
        // 结束流体求解end solveFlow
        Info << "\nmixture.correct" << endl;
        mixture.correct();
        turbulence->correct();

        // #include "alphaCourantNo.H"
        // #include "setDeltaT.H"
        // mixture.correct();
        
        
        Info << "CalcVelocityCorrection" << endl;
        particleCloud.calcVelocityCorrection
        (
            rho,
            p,
            U,
            phiIB,
            voidfraction,
            udivmid,
            usdomain,
            cella,
            alpha1
        );
        p_rgh = p - rho*gh;
        // udivlater = fvc::div(U);

        // rho = rhof*voidfraction + rhop*(1 - voidfraction);
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << " ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        Info<< "--------------------------------------------------" << endl;
    }
    return 0;
}
// ************************************************************************* //