# dimension-inference-with-betti-curve-comparison
This collection of code on this page implements a _semi-generalized_ version of the data analysis approach (to infer the intrinsic dimension of data) in the final part of my PhD dissertation (explaining the approach and implenting it to infer the dimension of the olfactory space of lab mice). 

Specifically, in my dissertation, the direction controls are **positive** and **uniform**. In addition to positive and uniform control, the code here also allows the direction to be sampled from a cone with a specified _cone angle_ &theta;.

## Prerequisites
Install the Simplicial package. Details can be found [here](https://github.com/nebneuron/Simplicial.jl).
```Julia
import Pkg; Pkg.add("Simplicial");
```
Install the JLD package. Details can be dound [here](https://github.com/JuliaIO/JLD.jl).
```Julia
import Pkg; Pkg.add("JLD");
```

## Basic Usage
Prepare your data matrix and save it as a .jld file; below, `M` is an m-by-n data matrix, where m is the number of sensors and n is sample size (i.e. the number of sample points). Make sure `M` is saved in the key `"response"`.
```Julia
using JLD
save("your_data_path.jld", "response", M)
```
Load the core code.
```Julia
using Simplicial
include("src/dim-inf-with-betti-curve-comparison.jl")
```
Modify the parameters in `run_unified_real_fvec_PI_bettis.jl` and `run_unified_ctrl_fvec_PI_bettis_parallel.jl`, and run the code below.
```Julia
include("run_unified_real_fvec_PI_bettis.jl")
include("run_unified_ctrl_fvec_PI_bettis_parallel.jl")
```
