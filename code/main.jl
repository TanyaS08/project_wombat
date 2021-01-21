## This will activate the environment in the code sub-folder
import Pkg
Pkg.activate("code")

## We then load the packages
using EcologicalNetworks

## Import the functions and methods we need
include(joinpath(pwd(), "code", "lib", "main.jl"))
