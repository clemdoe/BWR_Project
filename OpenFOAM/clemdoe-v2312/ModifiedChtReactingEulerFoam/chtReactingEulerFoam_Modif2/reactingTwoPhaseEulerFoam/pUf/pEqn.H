surfaceScalarField alphaf1("alphaf1", fvc::interpolate(alpha1));
surfaceScalarField alphaf2("alphaf2", scalar(1) - alphaf1);

surfaceScalarField alphaRhof10
(
    "alphaRhof10",
    fvc::interpolate
    (
        max(alpha1.oldTime(), phase1.residualAlpha())
       *rho1.oldTime()
    )
);

surfaceScalarField alphaRhof20
(
    "alphaRhof20",
    fvc::interpolate
    (
        max(alpha2.oldTime(), phase2.residualAlpha())
       *rho2.oldTime()
    )
);

// Drag coefficient
const surfaceScalarField Kdf("Kdf", fluid.Kdf());

// Diagonal coefficients
PtrList<surfaceScalarField> AFfs(fluid.AFfs());

PtrList<surfaceScalarField> rAUfs;
rAUfs.append
(
    new surfaceScalarField
    (
        IOobject::groupName("rAUf", phase1.name()),
        1.0
       /(
            byDt(alphaRhof10)
          + fvc::interpolate(U1Eqn.A())
          + AFfs[0]
        )
    )
);
rAUfs.append
(
    new surfaceScalarField
    (
        IOobject::groupName("rAUf", phase2.name()),
        1.0
       /(
            byDt(alphaRhof20)
          + fvc::interpolate(U2Eqn.A())
          + AFfs[0]
        )
    )
);
const surfaceScalarField& rAUf1 = rAUfs[0];
const surfaceScalarField& rAUf2 = rAUfs[1];

AFfs.clear();

// Explicit force fluxes
PtrList<surfaceScalarField> phiFfs(fluid.phiFfs(rAUfs));
const surfaceScalarField& phiFf1 = phiFfs[0];
const surfaceScalarField& phiFf2 = phiFfs[1];

// OJOOOO!!!!! EStoo es una chanchada!! pero lo estoy probando proque se me trababa
while (pimple.correct())
//for(int i=0; i<=2;i++)
{
    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    surfaceScalarField rhof1(fvc::interpolate(rho1));
    surfaceScalarField rhof2(fvc::interpolate(rho2));

    // Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF.correctBoundaryFlux(U1, phi1);
    MRF.correctBoundaryFlux(U2, phi2);

    const surfaceScalarField alpharAUf1
    (
        IOobject::groupName("alpharAUf", phase1.name()),
        max(alphaf1, phase1.residualAlpha())*rAUf1
    );

    const surfaceScalarField alpharAUf2
    (
        IOobject::groupName("alpharAUf", phase2.name()),
        max(alphaf2, phase2.residualAlpha())*rAUf2
    );

    // Combined buoyancy and force fluxes
    const surfaceScalarField ghSnGradRho
    (
        "ghSnGradRho",
        ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    const surfaceScalarField phigF1
    (
        alpharAUf1
       *(
            ghSnGradRho
          - alphaf2*(rhof1 - rhof2)*(g & mesh.Sf())
        )
      + phiFf1
    );

    const surfaceScalarField phigF2
    (
        alpharAUf2
       *(
            ghSnGradRho
          - alphaf1*(rhof2 - rhof1)*(g & mesh.Sf())
        )
      + phiFf2
    );

    // Predicted fluxes
    surfaceScalarField phiHbyA1
    (
        IOobject::groupName("phiHbyA", phase1.name()),
        phi1
    );

    phiHbyA1 =
        rAUf1
       *(
            alphaRhof10*byDt(MRF.absolute(phi1.oldTime()))
          + fvc::flux(U1Eqn.H())
        )
      - phigF1;

    surfaceScalarField phiHbyA2
    (
        IOobject::groupName("phiHbyA", phase2.name()),
        phi2
    );

    phiHbyA2 =
        rAUf2
       *(
            alphaRhof20*byDt(MRF.absolute(phi2.oldTime()))
          + fvc::flux(U2Eqn.H())
        )
      - phigF2;

    // Drag fluxes
    PtrList<surfaceScalarField> phiKdPhifs(fluid.phiKdPhifs(rAUfs));
    const surfaceScalarField& phiKdPhif1 = phiKdPhifs[0];
    const surfaceScalarField& phiKdPhif2 = phiKdPhifs[1];

    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf1*(phiHbyA1 - phiKdPhif1) + alphaf2*(phiHbyA2 - phiKdPhif2)
    );

    MRF.makeRelative(phiHbyA);

    phiKdPhifs.clear();

    // Construct pressure "diffusivity"
    const surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf1*alpharAUf1 + alphaf2*alpharAUf2)
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf1.boundaryField()*phi1.boundaryField()
              + alphaf2.boundaryField()*phi2.boundaryField()
            )
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp1, pEqnComp2;
    if (phase1.compressible())
    {
        if (pimple.transonic())
        {
            surfaceScalarField phid1
            (
                IOobject::groupName("phid", phase1.name()),
                fvc::interpolate(psi1)*phi1
            );

            pEqnComp1 =
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + correction
                (
                    (alpha1/rho1)*
                    (
                        psi1*fvm::ddt(p_rgh)
                      + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
                    )
                );

            pEqnComp1.ref().relax();
        }
        else
        {
            pEqnComp1 =
                (
                    fvc::ddt(alpha1, rho1) + fvc::div(phase1.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), rho1)
                )/rho1
              + (alpha1*psi1/rho1)*correction(fvm::ddt(p_rgh));
        }
    }
    if (phase2.compressible())
    {
        if (pimple.transonic())
        {
            surfaceScalarField phid2
            (
                IOobject::groupName("phid", phase2.name()),
                fvc::interpolate(psi2)*phi2
            );

            pEqnComp2 =
                (
                    fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + correction
                (
                    (alpha2/rho2)*
                    (
                        psi2*fvm::ddt(p_rgh)
                      + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
                    )
                );

            pEqnComp2.ref().relax();
        }
        else
        {
            pEqnComp2 =
                (
                    fvc::ddt(alpha2, rho2) + fvc::div(phase2.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), rho2)
                )/rho2
              + (alpha2*psi2/rho2)*correction(fvm::ddt(p_rgh));
        }
    }

    // Add option sources
    {
        if (fvOptions.appliesToField(rho1.name()))
        {
            tmp<fvScalarMatrix> optEqn1 = fvOptions(alpha1, rho1);
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= (optEqn1 & rho1)/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- (optEqn1 & rho1)/rho1, p_rgh);
            }
        }
        if (fvOptions.appliesToField(rho2.name()))
        {
            tmp<fvScalarMatrix> optEqn2 = fvOptions(alpha2, rho2);
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= (optEqn2 & rho2)/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- (optEqn2 & rho2)/rho2, p_rgh);
            }
        }
    }

    // Add mass transfer
    {
        PtrList<volScalarField> dmdts(fluid.dmdts());
        if (dmdts.set(0))
        {
            if (pEqnComp1.valid())
            {
                pEqnComp1.ref() -= dmdts[0]/rho1;
            }
            else
            {
                pEqnComp1 = fvm::Su(- dmdts[0]/rho1, p_rgh);
            }
        }
        if (dmdts.set(1))
        {
            if (pEqnComp2.valid())
            {
                pEqnComp2.ref() -= dmdts[1]/rho2;
            }
            else
            {
                pEqnComp2 = fvm::Su(- dmdts[1]/rho2, p_rgh);
            }
        }
    }

    // Cache p prior to solve for density update
    volScalarField p_rgh_0("p_rgh_0", p_rgh);

    //OJOO!!!! Esto probando esto pero es una chanchada
    while (pimple.correctNonOrthogonal())
    //for(int i=0; i<=3;i++)
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        {
            fvScalarMatrix pEqn(pEqnIncomp);

            if (pEqnComp1.valid())
            {
                pEqn += pEqnComp1();
            }

            if (pEqnComp2.valid())
            {
                pEqn += pEqnComp2();
            }

            pEqn.solve();
        }
        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        //if (i==3)
        {
            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            phi = phiHbyA + pEqnIncomp.flux();

            const surfaceScalarField phi1s
            (
                phiHbyA1 + alpharAUf1*mSfGradp
            );

            const surfaceScalarField phi2s
            (
                phiHbyA2 + alpharAUf2*mSfGradp
            );

            surfaceScalarField phir
            (
                ((phi2s + rAUf2*Kdf*phi1s) - (phi1s + rAUf1*Kdf*phi2s))
               /(1 - rAUf1*rAUf2*sqr(Kdf))
            );

            phi1 = phi - alphaf2*phir;
            phi2 = phi + alphaf1*phir;

            U1 = fvc::reconstruct(MRF.absolute(phi1));
            U1.correctBoundaryConditions();
            fvOptions.correct(U1);

            U2 = fvc::reconstruct(MRF.absolute(phi2));
            U2.correctBoundaryConditions();
            fvOptions.correct(U2);

            // Set the phase dilatation rates
            if (pEqnComp1.valid())
            {
                phase1.divU(-pEqnComp1 & p_rgh);
            }
            if (pEqnComp2.valid())
            {
                phase2.divU(-pEqnComp2 & p_rgh);
            }
        }
    }

    // Update and limit the static pressure
    p = max(p_rgh + rho*gh, pMin);

    // Limit p_rgh
    p_rgh = p - rho*gh;

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh - p_rgh_0);
    rho2 += psi2*(p_rgh - p_rgh_0);

    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
}
