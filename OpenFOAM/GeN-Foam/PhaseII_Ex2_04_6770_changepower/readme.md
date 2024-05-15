Phase II Exercise 2 is in the steady-state DNB condition. In these cases, the power level increases gradually to the critical heat flux condition.

In this exercise, we use similar boundary conditions (defined in 0/fluidRegion folder), the same model (defined in constant/fluidRegion folder), and the same solution strategy (defined in system folder).

This exercise used the distribution both in axial and radial directions. So we added a powerDensity.fixedPower file in 0/fluidRegion folder to represent the initial axial power distribution. Then we defined the power increase ratio against the initial value in each time position (See the file in constant/fluidRegion/phaseProperties - powerTimeProfile block).
In this exercise, we need to calculate the process of single-phase liquid to two-phase water, and then to critical heat flux condition. So we implemented a CHF look-up table model to predict the CHF value.

We can get the results of DNB power and DNB axial position in each case and compare them with experimental values.