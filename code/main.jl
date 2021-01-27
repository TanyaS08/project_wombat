## This will activate the environment in the code sub-folder
import Pkg
Pkg.activate("code")

## We then load the packages
using EcologicalNetworks

## Import the functions and methods we need
include(joinpath(pwd(), "code", "lib", "main.jl"))

#=

NOTE refer back to Fortin & Dale (2005) and Barbujani (1989) when you inevitibly get stuck

=#

## STEP 1: Create the network surface

#=

This means theoretically that we would need the values per a
network/community so that they can be situated in space (co-ords) as
well as an actual measure of the community

=#

## STEP 2: Connecting the networks in space

# assuming points are in a lattice configuration (uniformly contiguous)

# normalise co-ords so that they are between 0 - 1

# See Delaunay Triangulation (I have no idea why I said this...)

## STEP 3: Traverse and 'collate' the surface

#= RATE OF CHANGE: first partial derivative of a variable in
 𝑋 and 𝑌 spatial direction among four (4) adjacent samples/squares

        𝑚 = √[∂𝑓(𝑋,𝑌)/∂𝑋]²+ [∂𝑓(𝑋,𝑌)/∂𝑌]²

        where 𝑓(𝑋,𝑌) is a bilinear function of values 𝑧ᵢ at the
        four sampling locations (𝑖 = 1, 2, 3, 4)
        here 𝑧ᵢ would be our vaariable of interest e.g. a network property

        𝑓(𝑋,𝑌) = 𝑧₁(1 - 𝑋)(1 - 𝑌) + 𝑧₂𝑋(1 - 𝑌) + 𝑧₃𝑋𝑌 + 𝑧₄(1 - 𝑋)𝑌

        and

        ∂𝑓(𝑋,𝑌)/∂𝑋 = 𝑧₂ - 𝑧₁ + 𝑋(𝑧₁ - 𝑧₂ + 𝑧₃ - 𝑧₄)

        and

        ∂𝑓(𝑋,𝑌)/∂𝑌 = 𝑧₃ - 𝑧₁ + 𝑋(𝑧₁ - 𝑧₂ + 𝑧₃ - 𝑧₄)

    𝑚 it the 'new centroid' and represents the 'average' of 4 samples/squares

=#

## STEP 4: 'Quantify' the gradients

#= i.e. the rate of chnage so orientation and angle of the CHANGE

    θ = tan⁻¹[(∂𝑓/∂𝑥)/(∂𝑓/∂𝑦)] + Δ

    where

    Δ = 0°, if (∂𝑓/∂𝑥) =/> 0
        OR 180°, otherwise

=#

## STEP 5: Threshold values

#= i.e. deciding what would define a 'candidate' boundary

    - Rank rates of change descending order
       - NOTE same rate should be same percintile class even if at cutoff
    - Split into categories/groups based on chosen threshold and candidate
      boundaries
        e.g. if using 10ᵗʰ percentile then top candidates would be the first 𝑋
             values where 𝑋 is 10% of the number of candidate boundaries
             Candidate boundaries = (𝑚 - 1)(𝑛 - 1)

    - Second order derivative (identify inflection point i.e. where boundary ends)
        i.e. locating boundaries within boundaries so
        Candidate boundaries = (𝑚 - 2)(𝑛 - 2)
=#

## STEP 6: Signif of candidiate boundaries

#= using a restricted randomized test (binomial test)

    - if using 10ᵗʰ percintile each boundary has 𝑝 = 0.1 of being a candidate boundary
    - if 𝑎 variables out of 𝑏 variables are candidate boundaries at a given locality
      then boundary signif

      Pr(𝑎|𝑏) = (𝑏 𝑎)0.1ᵃ0.9ᵇ⁻ᵃ

        where (𝑏 𝑎) is the number of ways to choose 𝑎 elements out of 𝑏

        if </= to 5% then considered significant

=#


## STEP 6: Test of boundaries are connected

#=
    See boundary statistics
=#

## Triangulation Wombling for non-lattice structured data

#=
    See Fortin & Dale (2005) when you think you're ready
=#

## STEP ?: Could we quantify 'multivariate' gradients... i.e. repeating
# the process for multiple properties of our networks and somehow
# comparing and contrasting
