for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    //Info << "------------Obtengo los coeficiente heat trasnfer para las ecuaciones de energia-----------" << endl;

    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();

    //Info << "------------Resuelvo la ecuacion de Energia del Vapor-----------" << endl;
    if (!phase1.isothermal())
    {
        fvScalarMatrix E1Eqn
        (
            phase1.heEqn()
         ==
           *heatTransfer[phase1.name()]
          + alpha1*rho1*(U1&g)
          + fvOptions(alpha1, rho1, thermo1.he())
        );

        E1Eqn.relax();
        fvOptions.constrain(E1Eqn);
        E1Eqn.solve();
        fvOptions.correct(thermo1.he());
    }


    //Info << "------------Resuelvo la ecuacion de Energia del Agua-----------" << endl;
    if (!phase2.isothermal())
    {
        fvScalarMatrix E2Eqn
        (
            phase2.heEqn()
         ==
           *heatTransfer[phase2.name()]
          + alpha2*rho2*(U2&g)
          + fvOptions(alpha2, rho2, phase2.thermoRef().he())
        );

        E2Eqn.relax();
        fvOptions.constrain(E2Eqn);
        E2Eqn.solve();
        fvOptions.correct(thermo2.he());
    }

    fluid.correctThermo();
}
