using SpatialEcology
using LinearAlgebra
using DataFrames
import CSV
using SimpleSDMLayers
using Delaunay
using StatsPlots
using Statistics

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

# Map
plot(
    convert(Float32, A);
    c=:batlow,
    clim=(0, maximum(A)),
    frame=:box,
    title="Species richness",
)

X = copy(A.grid)
replace!(X, nothing => 0)
X = convert(Matrix{Int64}, X)

# Matrices for the strength and gradient
𝑀 = zeros(Float64, size(X) .- 1)
Θ = copy(𝑀)

for j in 1:(size(X, 2) - 1), i in 1:(size(X, 1) - 1)
    𝑀[i, j], Θ[i, j] = _rateofchange(X[i:(i + 1), j:(j + 1)])
end

change = SimpleSDMResponse(𝑀, A)
for lat in latitudes(change)
    for lon in longitudes(change)
        if isnothing(A[lon, lat])
            change[lon, lat] = 0.0
        end
    end
end
replace!(change, 0.0 => nothing)

plot(change; frame=:box, title="Rate of change")

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
    plot!(Shape(x, y); lab="", lw=0, fill_z=log(m), aspectratio=1)
end
xaxis!("Longitude", extrema(longitudes(A)))
yaxis!("Latitude", extrema(latitudes(A)))
title!("Delaunay triangulation")

# Example with lattice and bioclim data
# This example is a little bit faster because it has a loop to avoid the squares with empty values
A = worldclim(1; left=-180.0, right=180.0, bottom=-62.0, top=90.0)
rescale!(A, (0.0, 1.0))

# Matrices for the strength and gradient
𝑀 = convert(Matrix{Union{Float32,Nothing}}, zeros(Float32, size(A) .- 1))
Θ = copy(𝑀)

for j in 1:size(𝑀, 2), i in 1:size(𝑀, 1)
    tmp = A.grid[i:(i + 1), j:(j + 1)]
    if !any(isnothing.(tmp))
        𝑀[i, j], Θ[i, j] = _rateofchange(convert(Matrix{eltype(A)}, tmp))
    else
        𝑀[i, j], Θ[i, j] = (nothing, nothing)
    end
end

change = SimpleSDMResponse(𝑀, A)
angle = SimpleSDMResponse(Θ, A)

plot(change; dpi=600, c=:batlow)

# Tenth percentile but on the log of the rate of change
plot(rescale(log(change), [0.0, 0.90, 1.0]), dpi=600, c=:lapaz, legend=false)
title!("Possible boundaries")

qc = rescale(change, collect(0.0:0.01:1.0))

plot(
    change;
    frame=:box,
    title="Rate of change",
    clim=Tuple(quantile(collect(change), [0.1, 0.90])),
    c=:roma,
    dpi=600,
)

plot(qc; c=:roma, dpi=600, title="Quantiles of the rate of change")