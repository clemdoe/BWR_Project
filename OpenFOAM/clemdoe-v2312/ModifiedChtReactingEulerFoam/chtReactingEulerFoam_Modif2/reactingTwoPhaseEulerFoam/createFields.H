#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<twoPhaseSystem> fluidPtr
(
    twoPhaseSystem::New(mesh)
);
twoPhaseSystem& fluid = fluidPtr();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = fluid.phase1().thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());



// Esto lo agregue yo para computar el flujo de agua y vapor en particular a la salida
Info<< "Cargando el campo phiSteam\n" << endl;

    surfaceScalarField phiSteam
    (
        IOobject
        (
            "phiSteam",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(fluid.phase1().thermoRef().rho()*fluid.phase1().URef())
            & mesh.Sf()
    );

Info<< "Cargando el campo phiWater\n" << endl;

surfaceScalarField phiWater
(
    IOobject
    (
        "phiWater",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(fluid.phase2().thermoRef().rho()*fluid.phase2().URef())
        & mesh.Sf()
);
