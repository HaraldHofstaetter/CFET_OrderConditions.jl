# CFET_OrderConditions.jl

## Requirements
`CFET_OrderConditions.jl` requires the package `Giac`, a Julia interface to the Giac computer algebra system, see https://github.com/HaraldHofstaetter/Giac.jl.


## Installation
```julia
julia> Pkg.clone("https://github.com/HaraldHofstaetter/CFET_OrderConditions.jl")
```

## Examples
To get easy access to the examples, copy them into the home directory:
```julia
julia> cp(joinpath(homedir(), ".julia/v0.4/CFET_OrderConditions/examples/"), joinpath(homedir(), "CFET_OrderCondition_examples"), remove_destination=true)
```
Then 'Expokit_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [CFET_OrderConditions_examples.ipynb](https://github.com/HaraldHofstaetter/CFET_OrderConditions.jl/blob/master/examples/CFET_OrderConditions.ipynb)
