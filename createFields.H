Info<< "========Createfields for myinterid 2.0========" << endl;
#include "createRDeltaT.H"

Info<< "Reading voidfraction" << endl;
volScalarField voidfraction
(
    IOobject
    (
        "voidfraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
voidfraction.oldTime();

volScalarField alphamax
(
    IOobject
    (
        "alphamax",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)
);

volScalarField alphamin
(
    IOobject
    (
        "alphamin",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
);

// voidfraction.oldTime();

Info<< "Reading field p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Info<< "Reading field p_rghpre" << endl;
// volScalarField p_rghpre
// (
//     IOobject
//     (
//         "p_rghpre",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     p_rgh
// );

Info<< "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Create face flux field phi" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

// volScalarField divU
// (
//     IOobject
//     (
//         "divU",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     fvc::div(phi)
// );

Info<< "Reading particle velocity field Us" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties" << endl;

// 创建mixture对象，mixture中包含了"alpha1"，"alpha2"两个体心标量场
// multiphaseMixture mixture(U, phi);
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
    // voidfraction*rhof + (1 - voidfraction)*rhop
);

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// - overall non-temporal flux
surfaceScalarField phiAlpha
(
    IOobject
    (
        "phiAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U*alpha1) & mesh.Sf()
);

// Construct incompressible turbulence model

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<< "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PISO"),
    pRefCell,
    pRefValue
);

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

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"
// ==========================================================================
// variables for particleCloud.evolve()
Info<< "Reading field interFace" << endl;
volScalarField interFace
(
    IOobject
    (
        "interFace",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
);
Info<< "Reading field ibdrag" << endl;
volVectorField ibdrag
(
    IOobject
    (
        "ibdrag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ibd", dimensionSet(1,-2,-2,0,0,0,0),vector(1.0,0.0,0.0))
);
volVectorField ibdragP
(
    IOobject
    (
        "ibdragP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ibdp", dimensionSet(1,-2,-2,0,0,0,0),vector(1.0,0.0,0.0))
);
volVectorField ibdragV
(
    IOobject
    (
        "ibdragV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ibdv", dimensionSet(1,-2,-2,0,0,0,0),vector(1.0,0.0,0.0))
);
// =============================================================================
// variables for particleCloud.calcVelocityCorrection()
Info<< "Reading field phiIB" << endl;
volScalarField phiIB
(
    IOobject
    (
        "phiIB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField udivmid
(
    IOobject
    (
        "udivmid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
);
volVectorField usdomain
(
    IOobject
    (
        "usdomain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector(0,0,0))
);

volScalarField cella
(
    IOobject
    (
        "cella",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
);
// ============================================================================

IOobject alphaPhi10Header
(
    "alphaPhi10",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);
// const bool alphaRestart = alphaPhi10Header.typeHeaderOk<surfaceScalarField>();
// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    IOobject
    (
        "alphaPhi10",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
//            IOobject::AUTO_WRITE
    ),
    linearInterpolate(U*alpha1) & mesh.Sf()
);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

Info<< "==========End createfields========\n" << endl;
