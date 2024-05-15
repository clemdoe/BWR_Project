# 1D PSBT SC

The OECD/NRC PWR PSBT benchmark was organized based on the NUPEC database. It is a well-known benchmark for code validation.


## Phase I Exercise 1 12223

This tutorial presents a simple 1-D case for water boiling based on the *OECD/NRC Benchmark based on NUPEC Pressurised Water Reactor (PWR) Subchannel and Bundle Tests (PSBT)*, case 12223 (exercise 1).

Please note that models for water boiling are still preliminary, incomplete (missing models for boiling crisis) and in Beta testing.


## Phase II

This section presents 3 exercises of *Pressurised Water Reactor (PWR) Subchannel and Bundle Tests (PSBT) - Phase II*.


### Exercise 1 01-5215

Phase II Exercise 1 is in the water boiling condition, which aims at calculating the liquid temperature at the specific position.

In this exercise, we use similar boundary conditions (defined in `0/fluidRegion` folder), the same model (defined in `constant/fluidRegion` folder), and the same solution strategy (defined in `system` folder). 

This exercise used a challenging power gradient in the radial direction. The model aimed at calculating the liquid temperature at the axial height of 4.775 m to compare with experimental values.


### Exercise 2 04-6770

Phase II Exercise 2 is in the steady-state DNB condition. In these cases, the power level increases gradually to the critical heat flux condition. 

In this exercise, we use similar boundary conditions (defined in `0/fluidRegion` folder), the same model (defined in `constant/fluidRegion` folder), and the same solution strategy (defined in `system` folder). 

This exercise used the distribution both in axial and radial directions. So we added a powerDensity.fixedPower file in `0/fluidRegion` folder to represent the initial axial power distribution. Then we defined the power increase ratio against the initial value in each time position (See the file in constant/fluidRegion/phaseProperties - powerTimeProfile block).

In this exercise, we need to calculate the process of single-phase liquid to two-phase water, and then to critical heat flux condition. So we implemented a CHF look-up table model to predict the CHF value.

We can get the results of DNB power and DNB axial position in each case and compare them with experimental values.


### Exercise 3 11-0312

Phase II Exercise 3 is in transient DNB condition. In these cases, the boundary condition is changing over time. Each case also gradually reaches the CHF condition through the change of boundary condition.

In this exercise, we use the different boundary conditions (defined in `0/fluidRegion` folder), the same model (defined in `constant/fluidRegion` folder), and the same solution strategy (defined in `system` folder). 

This exercise used the distribution both in axial and radial directions. So we added a powerDensity.fixedPower file in `0/fluidRegion` folder to represent the initial axial power distribution. Then we defined the power increase ratio against the initial value in each time position (See the file in `constant/fluidRegion/phaseProperties` - powerTimeProfile block).

This exercise is transition DNB cases. The boundary condition is changing over time.

In this exercise, we also need to calculate the process of single-phase liquid to two-phase water, and then the critical heat flux condition. So the CHF look-up table is also implemented in this exercise.

We can get the results of DNB power and the time of occurring DNB in each case and compare them with experimental values.


## Run simulation

All cases are provided with the same running script. Simply execute the following commands in the terminal.

```bash
./Allclean
./Allrun
# or
./Alltest
```
