//regionProperties rp(runTime);

const wordList scalarNames(rp["scalar"]);

PtrList<fvMesh> scalarRegions(scalarNames.size());

forAll(scalarNames, i)
{
    Info<<"------------ Nueva Region Agregada ------------ " << endl;
    Info<< "Create scalar mesh for region " << scalarNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    scalarRegions.set
    (
        i,
        new fvMesh
        (
            IOobject
            (
                scalarNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
