{
    volScalarField& he = thermo.he();

    //Info << "El campo h fluido inicial - max" << max(he.primitiveField()) << " - min - " << min(he.primitiveField()) <<endl;


    fvScalarMatrix EEqn
    (
        fvm::ddt(rhofl, he) + fvm::div(phifl, he)
      + fvc::ddt(rhofl, Kfl) + fvc::div(phifl, Kfl)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phifl/fvc::interpolate(rhofl), Ufl),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulencefl.alphaEff(), he)
     ==
        rhofl*(Ufl&gfl)
      + fvOptions(rhofl, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    //Info << "El campo h fluido final - max" << max(he.primitiveField()) << " - min - " << min(he.primitiveField()) <<endl;

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
