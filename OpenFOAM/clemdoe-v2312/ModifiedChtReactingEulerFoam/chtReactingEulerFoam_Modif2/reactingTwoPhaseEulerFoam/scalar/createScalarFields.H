// ---------------------------------------------------------------- //

// Initialise fluid field pointer lists
PtrList<volScalarField> Tscalar(scalarRegions.size());
PtrList<volScalarField> rhoScalar(scalarRegions.size());
PtrList<volVectorField> Uscalar(scalarRegions.size());
PtrList<surfaceScalarField> phiScalar(scalarRegions.size());
PtrList<fv::options> scalarFvOptions(scalarRegions.size());
PtrList<dimensionedScalar> DTscalar(scalarRegions.size());
PtrList<solidThermo> thermoScalar(scalarRegions.size());


PtrList<volScalarField> KScalar(scalarRegions.size());
//PtrList<volScalarField> dpdtFluid(scalarRegions.size());
PtrList<uniformDimensionedVectorField> gScalar(scalarRegions.size());


// Populate fluid field pointer lists
forAll(scalarRegions, i)
{
    Info<< "*** Reading scalar mesh thermophysical properties for region "
        << scalarRegions[i].name() << nl << endl;

    Info<< "    Adding to thermos\n" << endl;
    thermoScalar.set(i, solidThermo::New(scalarRegions[i]));


Info<< "    Adding to rhoFluid\n" << endl;
    rhoScalar.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                scalarRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoScalar[i].rho()
        )
    );


Info<< "    Adding to field Uscalar\n" << endl;

Uscalar.set
(
    i,
    new volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            scalarRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        scalarRegions[i]
    )
);


Info<< "    Adding/calculating face flux field phi\n" << endl;

phiScalar.set
(
    i,
    new surfaceScalarField
    (
    IOobject
    (
        "phi",
        runTime.timeName(),
        scalarRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(rhoScalar[i]*Uscalar[i])
           &  scalarRegions[i].Sf()

    )
);


Info<< "    Adding to gScalar\n" << endl;
gScalar.set
(
    i,
    new uniformDimensionedVectorField
    (
        IOobject
        (
            "g",
            runTime.constant(),
            scalarRegions[i],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
);

Info<< "    Adding to KScalar\n" << endl;
KScalar.set
(
    i,
    new volScalarField
    (
        "K",
        0.5*magSqr(Uscalar[i])
    )
);

/*Info<< "    Adding to dpdtScalar\n" << endl;
dpdtScalar.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            scalarRegions[i],
        ),
        scalarRegions[i],,
        dimensionedScalar
        (
            thermoFluid[i].p().dimensions()/dimTime,
            0
        )
    )
);*/


Info<< "    Adding fvOptions\n" << endl;
scalarFvOptions.set
(
    i,
    new fv::options(scalarRegions[i])
);

}
