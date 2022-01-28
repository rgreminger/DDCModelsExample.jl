# DDCModelsExample.jl

This package implements Julia code for the basic model in "Dynamic Discrete Choice Models: Methods, Matlab Code, and Exercises" of Abbring and Klein (2020). Detailed explanations of the model and the Matlab code are available here: [https://jabbring.github.io/dynamic-discrete-choice/dynamicDiscreteChoice.m.html](https://jabbring.github.io/dynamic-discrete-choice/dynamicDiscreteChoice.m.html)

The goal is to set up the basic model and offer a package structure that easily extends to more complex DDC problems. Most functions are explained through docstrings (e.g. type `?estimate_model` to see the documentation). Steps in individual functions are explained in the respective source files. 

The package implements the type `DDCbasic`; a structure for the basic model containing the parameters and simulation/estimation options. To check which functions and types the package exports, check the file *src/DDCModelsExample.jl*. 

Various functions implement multi-threading for faster execution times. To use e.g. 8 threads, start Julia with `julia --threads 8`. In VScode, the number of threads can be set in the Julia extension settings. Note that for a running Julia instance, the number of threads cannot be changed. 

To see how to use the package, check the file *test/runtests.jl*, which sets up unit tests based on Julia's `Test` framework. Specifically, it: 

1. constructs a `DDCbasic` model
2. simulates data from this model 
3. estimates the transition matrix 
4. estimates parameters 
5. calculates standard errors of the parameter estimates. 

## Installing the package
Note, the package was tested and should be run with Julia version 1.6. or newer. 
1. Download and extract the package from Github (e.g. through `git clone`)
2. In Julia, navigate to the package folder (or open the folder in VSCode)
3. Activate environment: `using Pkg; Pkg.activate(".")`
4. Install dependencies: `Pkg.instantiate()`
5. Test to verify it runs correctly: `Pkg.test("DDCModelsExample")`
6. Load the package with `using DDCModelsExample` 


## Files in *src*: 
- *types.jl*: Contains all type definitions 
- *support.jl*: Defines support functions that are useful for (potentially) various differing models 
- *models/DDCbasic.jl*: Functions specific to the `DDCbasic` model (e.g. flowpayoffs, loglikelihood etc.)
- *test/runtests.jl*: Implement unit tests
- *test/speedtest.jl*: Setup for running speed tests



