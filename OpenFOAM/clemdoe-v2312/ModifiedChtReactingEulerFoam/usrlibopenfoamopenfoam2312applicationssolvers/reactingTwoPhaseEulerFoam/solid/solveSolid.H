{
    while (pimpleSolid.correctNonOrthogonal())
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rhosol, hsol)
          - (
                thermosol.isotropic()
              ? fvm::laplacian(betav*thermosol.alpha(), hsol, "laplacian(alpha,h)")
              : fvm::laplacian(betav*taniAlpha(), hsol, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rhosol, hsol)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve();

        fvOptions.correct(hsol);
    }
}

thermosol.correct();

Info<< "Min/max T:" << min(thermosol.T()).value() << ' '
    << max(thermosol.T()).value()<<'\n' << endl;
