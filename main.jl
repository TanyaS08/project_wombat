# # README

# ![](assets/header.png){#fig:head}

import Pkg
Pkg.activate("/Users/tanyastrydom/Documents/Uni/project_wombat")


# this can be pruned at some point...

using Base
using CSV: CSV
using DataFrames
using Delaunay
using LinearAlgebra
using SimpleSDMLayers
using SimpleSDMLayers: convert
using SpatialEcology
using Statistics
using StatsBase
using StatsPlots
using Plots

theme(:mute)
default(; frame=:box)

# Import the functions and methods we need

include(joinpath(pwd(), "lib", "rateofchange.jl"))
include(joinpath(pwd(), "lib", "boundary.jl"))
include(joinpath(pwd(), "lib", "wombling.jl"))

# ## A Neutral landscapes example for lattice wombling

# First we create a 50 x 50 'cell' landscape. In this
# instance a planar landscape is it has a 'distinct'/sharp
# boundary 

using NeutralLandscapes
siz = 50, 50
A = Matrix(rand(EdgeGradient(), siz))

# Plot the landscape so we have an _a priori_ idea
# of what if looks like 

plot(heatmap(A))

savefig("assets/PlanarLandscape.png")

# ![](assets/PlanarLandscape.png){#fig:landscape}

# Calculate the rate of change (𝑀) and direction of 
# change (θ) using `wombling()`


womble = wombling(A)

plot(
    plot(heatmap(womble.M), title = "Rate of Change"),
    plot(heatmap(womble.Θ), title = "Direction of Change")
    )

savefig("assets/PlanarWomble.png")

# ![](assets/PlanarWomble.png){#fig:womble}

# We can extract the candidate boundaries using the 
# rate of change caluclated from `wombling()`

boundary = boundaries(womble.M)

plot(
    heatmap(boundary),
    title = "Candidate boundaries"
)

savefig("assets/PlanarBoundaries.png")

# ![](assets/PlanarBoundaries.png){#fig:womble}