    scalar CoNumFluid = -great;

    forAll(fluidRegions, regionI)
    {
        CoNumFluid = max
        (
            compressibleCourantNoFluid
            (
                fluidRegions[regionI],
                runTime,
                rhoFluid[regionI],
                phiFluid[regionI]
            ),
            CoNumFluid
        );
    }
