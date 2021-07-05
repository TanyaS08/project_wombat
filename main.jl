using SimpleSDMLayers: hcat
import Pkg
Pkg.activate("/Users/tanyastrydom/Documents/Uni/project_wombat")

using CSV: CSV
using DataFrames
using Delaunay
using LinearAlgebra
using NeutralLandscapes
using SimpleSDMLayers
using SpatialEcology
using Statistics
using StatsPlots
using Plots

theme(:mute)
default(; frame=:box)

#=
NOTE refer back to Fortin & Dale (2005) and Barbujani (1989) when you inevitably get stuck
=#

# Import the functions and methods we need
include(joinpath(pwd(), "lib", "rateofchange.jl"))

# Transform the amphibian data into a raster of richness
amphibians = DataFrame(
    CSV.File(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv"))
)
select!(amphibians, Not(:coords))

# Get the richness map by grouping on Lat and Lon, and summing the content of the cells
n_species = combine(
    groupby(stack(amphibians, Not([:Long, :Lat])), [:Long, :Lat]),
    :value => sum => :richness,
)

# Prepare the amphibian data with correct column names
rename!(n_species, :Long => :longitude, :Lat => :latitude)
sort!(n_species, [:latitude, :longitude])

lats = sort(unique(n_species.latitude))
lons = sort(unique(n_species.longitude))
Δlat = lats[2] - lats[1]
Δlon = lons[2] - lons[1]
G = fill(nothing, length(lats), length(lons))
G = convert(Matrix{Union{Nothing,Int64}}, G)
A = SimpleSDMResponse(
    G,
    minimum(lons) - 0.5Δlon,
    maximum(lons) + 0.5Δlon,
    minimum(lats) - 0.5Δlat,
    maximum(lats) + 0.5Δlat,
)


for site in eachrow(n_species)
    A[site.longitude, site.latitude] = site.richness
end

A = convert(Float64, A)

# Map
plot(
    convert(Float32, A);
    c=:batlow,
    clim=(0, maximum(A)),
    frame=:box,
    title="Species richness",
)

# Do a Delaunay thingie from the sites
amph_points = Matrix(n_species[!, [:longitude, :latitude]])
mesh = delaunay(amph_points)



plot()
for i in 1:size(mesh.simplices, 1)
    c = amph_points[mesh.simplices[i, :], :]
    x = c[:, 1]
    y = c[:, 2]
    z = n_species.richness[mesh.simplices[i, :]]
    m, t = _rateofchange(x, y, z)
    plot!(Shape(x, y); lab="", lw=0, fill_z=t, aspectratio=1)
end
xaxis!("Longitude", extrema(longitudes(A)))
yaxis!("Latitude", extrema(latitudes(A)))
title!("Delaunay triangulation")

c = hcat(x, y)[mesh.simplices[1, :], :]
x = c[:, 1]
y = c[:, 2]
z = n_species.richness[mesh.simplices[1, :]]
m, t = _rateofchange(x, y, z)



for i in 1:size(mesh.simplices, 1)
    c = amph_points[mesh.simplices[i, :], :]
    x = c[:, 1]
    y = c[:, 2]
    z = n_species.richness[mesh.simplices[i, :]]
    m, t = _rateofchange(x, y, z)
end

x = amph_points[:, 2]
y = amph_points[:, 1]
z = n_species.richness

hcat(x, y)

delaunay(hcat(x, y))

"""
    triangulationwombling(x::Vector{T}, y::Vector{T}, z::Vector) where {T<:Number}

Wrapper function that implements the triangualtion wombling algorithm for points 
that are irregualrly arranged in space.
"""
function triangulationwombling(x::Vector{T}, y::Vector{T}, z::Vector) where {T<:Number}
    length(x) >= 3  || throw(DimensionMismatch("x must have a minimum length of 3"))
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    length(x) == length(z) || throw(DimensionMismatch("x and z must have the same length"))

    # Do Delaunay thingie for sites
    mesh = delaunay(hcat(x, y))

    for i in 1:size(mesh.simplices, 1)
        c = hcat(x, y)[mesh.simplices[i, :], :]
        x = c[:, 1]
        y = c[:, 2]
        z = z[mesh.simplices[i, :]]
        return x, y, z
    end

    # Rate of change and direction
    return _rateofchange(x, y, z)
end

triangulationwombling(x, y, z)



# Example with lattice and bioclim data
# This example is a little bit faster because it has a loop to avoid the squares with empty values
A = worldclim(3; bottom=-60.0)
rescale!(A, (0.0, 1.0))
plot(A,
background_color = :transparent,
foreground_color = :black,)

png("figures/Isothermality")

# Matrices for the strength and gradient
𝑀 = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(A) .- 1))
Θ = copy(𝑀)

for j in 1:size(𝑀, 2), i in 1:size(𝑀, 1)
    tmp = A.grid[i:(i + 1), j:(j + 1)]
    if !any(isnothing.(tmp))
        tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀[i, j], Θ[i, j] = _rateofchange(tmp;
                                          X = stride(A, 1), Y = stride(A, 2)
                                          )
    else
        𝑀[i, j], Θ[i, j] = (nothing, nothing)
    end
end

change = SimpleSDMResponse(𝑀, A)
angle = SimpleSDMResponse(Θ, A)

# Colors for North, South, East, and West -- this is a square of complementary colors
CE, CS, CN, CW = colorant"#e3d96d", colorant"#714be3", colorant"#e35e40", colorant"#40e3a8"
plot(angle; 
    c=cgrad([CN, CE, CW, CS, CN], [0.0, 90.0, 180.0, 270.0, 360.0]), 
    clim=(0, 360.0), 
    dpi=400,
    background_color = :transparent,
    foreground_color = :black,)
title!("Direction of Change")

png("figures/DirectionofChange")

# Tenth percentile but on the log of the rate of change
plot(rescale(log(change), 
    [0.0, 0.90, 1.0]); 
    dpi=400, c=:lapaz, 
    legend=false,
    background_color = :transparent,
    foreground_color = :black,)
title!("Possible boundaries")

png("figures/PossibleBoundaries")

plot(
    log(change);
    frame=:box,
    title="Rate of change",
    dpi=400,
    background_color = :transparent,
    foreground_color = :black,
)

png("figures/RateofChange_log")

plot(
    rescale(change, collect(0.0:0.01:1.0));
    dpi=400,
    title="Quantiles of the rate of change",
    background_color = :transparent,
    foreground_color = :black,
)

png("figures/Quantiles_RateofChange")

# Example with lattice and NeutralLandscapes
# which is ugly because I didnt feel like thinking
siz = 50, 50


cluster = rand(RectangularCluster(4, 8), siz)
random = rand(NoGradient(), siz)
edge = rand(EdgeGradient(), siz)
planar = rand(PlanarGradient(), siz)
midpt = rand(MidpointDisplacement(0.75), siz)

A = cluster


𝑀 = convert(Matrix{Union{Float32}}, zeros(Float32, size(A) .- 1))
Θ = copy(𝑀)

for j in 1:size(𝑀, 2), i in 1:size(𝑀, 1)
    tmp = A[i:(i + 1), j:(j + 1)]
        #tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀[i, j], Θ[i, j] = _rateofchange(tmp)
end

𝑀_r = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(random) .- 1))
Θ_r = copy(𝑀_r)

for j in 1:size(𝑀_r, 2), i in 1:size(𝑀_r, 1)
    tmp = random[i:(i + 1), j:(j + 1)]
        #tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀_r[i, j], Θ_r[i, j] = _rateofchange(tmp)
end

𝑀_e = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(edge) .- 1))
Θ_e = copy(𝑀_e)

for j in 1:size(𝑀_e, 2), i in 1:size(𝑀_e, 1)
    tmp = edge[i:(i + 1), j:(j + 1)]
        #tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀_e[i, j], Θ_e[i, j] = _rateofchange(tmp)
end

𝑀_p = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(planar) .- 1))
Θ_p = copy(𝑀_p)

for j in 1:size(𝑀_p, 2), i in 1:size(𝑀_p, 1)
    tmp = planar[i:(i + 1), j:(j + 1)]
        #tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀_p[i, j], Θ_p[i, j] = _rateofchange(tmp)
end

𝑀_m = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(midpt) .- 1))
Θ_m = copy(𝑀_m)

for j in 1:size(𝑀_m, 2), i in 1:size(𝑀_m, 1)
    tmp = midpt[i:(i + 1), j:(j + 1)]
        #tmp = convert(Matrix{eltype(A)}, tmp)
        𝑀_m[i, j], Θ_m[i, j] = _rateofchange(tmp)
end

gr(color = :viridis, ticks = false, framestyle = :box, colorbar = false, aspectratio = 1)
plot(
    plot(heatmap(cluster), title = "Cluster"), heatmap(𝑀), heatmap(Θ),
    plot(heatmap(random), title = "Random"), heatmap(𝑀_r), heatmap(Θ_r),
    plot(heatmap(edge), title = "Edge"), heatmap(𝑀_e), heatmap(Θ_e),
    plot(heatmap(planar), title = "Planar"), heatmap(𝑀_p), heatmap(Θ_p),
    plot(heatmap(midpt), title = "Mid point"), heatmap(𝑀_m), heatmap(Θ_m),
    layout = (5, 3), size = (750, 1400))

png("figures/NeutralLandscapes")

plot(
    heatmap(Θ)
    )


sort(𝑀[:, 1], dims = 1)

C = Any[]

for i in 1:size(𝑀, 1)

    max = size(𝑀[:, i], 1)
    min = max - convert(Int64, round(max*0.1, digits = 0))
    a = 𝑀[:, i]
    b = partialsortperm(a, min:max)

    push!(C, a[b])

end

a = 𝑀[:, 1]
b = partialsortperm(a, 1:10)
collect(a[b])
