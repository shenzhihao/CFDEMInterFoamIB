{
    Info << "alphaEqu version 2.0" << endl;
    word alphaScheme("div(phi,alpha)");
    // Gauss vanLeer
    word alpharScheme("div(phirb,alpha)");
    // Gauss interfaceCompression
    // surfaceScalarField voidf = fvc::interpolate(voidfraction);
    // surfaceScalarField phivoid = phi * voidf;

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*mixture.nHatf());

    Info << "use voidfraction"<< endl;

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        // #include "alphaSuSp.H"
        // zeroField Su;
        // zeroField Sp;
        // zeroField divU;
        alphaPhi10 = 
        // fvc::flux(phivoid,alpha1,alphaScheme)
        fvc::flux(phi,alpha1,alphaScheme)
        // + fvc::interpolate(alpha1) * phivoidre
        // + fvc::flux(-fvc::flux(-phirvoid, alpha2, alpharScheme),alpha1,alpharScheme);
        + fvc::flux(-fvc::flux(-phir, alpha2, alpharScheme),alpha1,alpharScheme);
        
        // divU = fvc::div(phi);
        // divU = fvc::ddt(voidfraction);

        MULES::explicitSolve
        (
            geometricOneField(),
            // voidfraction,
            alpha1,
            phi,
            alphaPhi10,
            // Sp,
            zeroField(),
            // divU,
            // (Su + divU*min(alpha1(), scalar(1)))(),
            zeroField(),
            1,
            0
        );
        alpha1 = min(alpha1 , alphamax);
        alpha1 = max(alpha1 , alphamin);
        alpha2 = 1.0 - alpha1;
        mixture.correct();
        #include "rhofs.H"
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
    }
    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
