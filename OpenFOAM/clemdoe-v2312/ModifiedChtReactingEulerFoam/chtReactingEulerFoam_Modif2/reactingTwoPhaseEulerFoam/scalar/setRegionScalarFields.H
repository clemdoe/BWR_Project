    //const fvMesh& mesh = scalarRegions[i];

    surfaceScalarField& phiscal = phiScalar[i];

    volVectorField& Uscal = Uscalar[i];

    //Para usar esto tenes que cargar una thermofisical properties
    solidThermo& thermoScal = thermoScalar[i];

    volScalarField& hscal = thermoScal.he();

    volScalarField& rhoscal = rhoScalar[i];

    volScalarField& Kscal = KScalar[i];
    //volScalarField& dpdt = dpdtScalar[i];
    //volScalarField& p = thermo.p();
    //volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& gscal = gScalar[i];


    fv::options& fvOptions = scalarFvOptions[i];

    //pimpleNoLoopControl& pimple = pimples.pimple(i);

