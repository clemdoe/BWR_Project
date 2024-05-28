    const fvMesh& mesh = fluidRegions[i];

    rhoReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rhofl = rhoFluid[i];
    volVectorField& Ufl = UFluid[i];
    surfaceScalarField& phifl = phiFluid[i];

    compressible::turbulenceModel& turbulencefl = turbulenceFluid[i];
    volScalarField& Kfl = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& pfl = thermo.p();
    const volScalarField& psifl = thermo.psi();

    volScalarField& p_rghfl = p_rghFluid[i];

    const dimensionedVector& gfl = gFluid[i];
    const volScalarField& ghfl = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimpleFluid = pimples.pimple(i);

    pressureControl& pressureControlFluido = pressureControlFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];
