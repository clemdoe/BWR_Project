    // Resuelo la ecuacion de momento

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rhofl, Ufl) + fvm::div(phifl, Ufl)
      + turbulencefl.divDevRhoReff(Ufl)
     ==
        fvOptions(rhofl, Ufl)
    );
    fvVectorMatrix& UEqnFluido = tUEqn.ref();

    UEqnFluido.relax();

    fvOptions.constrain(UEqnFluido);

    if (pimpleFluid.momentumPredictor())
    {
        solve
        (
            UEqnFluido
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhofl)
                  - fvc::snGrad(p_rghfl)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(Ufl);
        Kfl = 0.5*magSqr(Ufl);
    }

    fvOptions.correct(Ufl);
