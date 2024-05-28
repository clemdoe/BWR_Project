
volScalarField rAUFluido("rAUFluido", 1.0/UEqnFluido.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rhofl*rAUFluido));
volVectorField HbyAFluido(constrainHbyA(rAUFluido*UEqnFluido.H(), Ufl, p_rghfl));

if (pimpleFluid.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rhofl)*mesh.magSf());

surfaceScalarField phiHbyAFluido
(
    "phiHbyAFluido",
    fvc::flux(rhofl*HbyAFluido)
);

const bool closedVolume = mesh.steady() && adjustPhi(phiHbyAFluido, Ufl, p_rghfl);
const bool adjustMass = closedVolume && !thermo.incompressible();

phiHbyAFluido += phig;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rghfl, rhofl, Ufl, phiHbyAFluido, rhorAUf, NullMRF());

{
    fvScalarMatrix p_rghEqnComp
    (
        fvc::ddt(rhofl) + psifl*correction(fvm::ddt(p_rghfl))
     ==
        fvOptions(psifl, p_rghfl, rhofl.name())
    );

    if (pimpleFluid.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psifl)/fvc::interpolate(rhofl))*phiHbyAFluido
        );

        phiHbyAFluido -= fvc::interpolate(psifl*p_rghfl)*phiHbyAFluido/fvc::interpolate(rhofl);

        p_rghEqnComp += fvm::div(phid, p_rghfl);
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    tmp<volScalarField> psip0(mesh.steady() ? tmp<volScalarField>() : psifl*pfl);

    while (pimpleFluid.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyAFluido)
          - fvm::laplacian(rhorAUf, p_rghfl)
        );

        fvScalarMatrix p_rghEqn(p_rghEqnComp + p_rghEqnIncomp);

        //Ver si despues lo agrego
        p_rghEqn.setReference
        (
            pressureControlFluido.refCell(),
            pressureControlFluido.refValue()
        );


        p_rghEqn.solve();

        if (pimpleFluid.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phifl = phiHbyAFluido + p_rghEqn.flux();


            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            Ufl = HbyAFluido
                + rAUFluido*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rhorAUf);
            Ufl.correctBoundaryConditions();
            fvOptions.correct(Ufl);
            Kfl = 0.5*magSqr(Ufl);
        }
    }

    pfl = p_rghfl + rhofl*ghfl;

    //Ver si despues lo agrego a esto
    //pressureControl.limit(pfl);

    // Thermodynamic density update
    if (!mesh.steady())
    {
        thermo.correctRho(psifl*pfl - psip0);
    }
}


// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(pfl);
}

//Esto lo hice porque rho se tiene q llamar asi
volScalarField rho = rhofl;


// Solve continuity
if (!mesh.steady())
{
    #include "rhoFluid.H"
    #include "compressibleContinuityErrs.H"
}

else
{
    #include "incompressible/continuityErrs.H"
}


// For closed-volume compressible cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass)
{
    pfl += (initialMass - fvc::domainIntegrate(rhofl))
        /fvc::domainIntegrate(psifl);
    p_rghfl = pfl - rhofl*ghfl;
    pfl.correctBoundaryConditions();
}

Info<< "Min/max rho:" << min(rhofl).value() << ' '
    << max(rhofl).value() << endl;



