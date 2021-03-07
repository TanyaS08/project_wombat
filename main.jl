using SpatialEcology
using LinearAlgebra
using CSV, DataFrames
using Statistics
using SimpleSDMLayers
using Delaunay
using StatsPlots

theme(:mute)
default(frame = :box)

#=
NOTE refer back to Fortin & Dale (2005) and Barbujani (1989) when you inevitably get stuck
=#

# Import the functions and methods we need
include(joinpath(pwd(), "lib", "rateofchange.jl"))


# Transform the amphibian data into a raster of richness
amphibians = DataFrame(CSV.File(joinpath(dirname(pathof(SpatialEcology)), "..", "data", "amph_Europe.csv")))
select!(amphibians, Not(:coords))
richness = combine(groupby(stack(amphibians, Not([:Long, :Lat])), [:Long, :Lat]), :value => sum => :richness)
rename!(richness, :Long => :longitude, :Lat => :latitude)
sort!(richness, [:latitude, :longitude])

lats = sort(unique(richness.latitude))
lons = sort(unique(richness.longitude))
Δlat = lats[2]-lats[1]
Δlon = lons[2]-lons[1]
G = fill(nothing, length(lats), length(lons))
G = convert(Matrix{Union{Nothing,Int64}}, G)
A = SimpleSDMResponse(G, minimum(lons)-0.5Δlon, maximum(lons)+0.5Δlon, minimum(lats)-0.5Δlat, maximum(lats)+0.5Δlat)

for site in eachrow(richness)
    A[site.longitude, site.latitude] = site.richness
end

# Map
plot(convert(Float32, A), c=:batlow, clim=(0, maximum(A)), frame=:box, title="Species richness")

X = copy(A.grid)
replace!(X, nothing => 0)
X = convert(Matrix{Int64}, X)

# Matrices for the strength and gradient
𝑀 = zeros(Float64, size(X).-1)
Θ = similar(𝑀)

for j in 1:(size(X,2)-1), i in 1:(size(X,1)-1)
    𝑀[i,j], Θ[i,j] = _rateofchange(X[i:(i+1),j:(j+1)])
end

change = SimpleSDMResponse(𝑀, A)
for lat in latitudes(change)
    for lon in longitudes(change)
        if isnothing(A[lon,lat])
            change[lon,lat] = 0.0
        end
    end
end
replace!(change, 0.0 => nothing)

plot(change, frame=:box, title="Rate of change")

# Do a Delaunay thingie from the sites
coords = Matrix(richness[!,[:longitude,:latitude]])
mesh = delaunay(coords)

plot()
for i in 1:size(mesh.simplices,1)
    c = coords[mesh.simplices[i,:],:]
    x = c[:,1]
    y = c[:,2]
    z = richness.richness[mesh.simplices[i,:]]
    m, t = _rateofchange(x, y, z)
    plot!(Shape(x,y), lab="", lw=0, fill_z = log(m), aspectratio=1)
end
xaxis!("Longitude", extrema(longitudes(A)))
yaxis!("Latitude", extrema(latitudes(A)))
title!("Delaunay triangulation")


# Example with lattice and bioclim data
A = worldclim(12; left=-180.0, right=180.0, bottom=-62.0, top=90.0)
rescale!(A, (0.0, 1.0))
plot(A)


X = copy(A.grid)
replace!(X, nothing => 0.0)
X = convert(Matrix{Float64}, X)

# Matrices for the strength and gradient
𝑀 = zeros(Float64, size(X).-1)
Θ = similar(𝑀)

for j in 1:(size(X,2)-1), i in 1:(size(X,1)-1)
    tmp = X[i:(i+1),j:(j+1)]
    if !any(isnothing.(tmp))
        if sum(tmp) != 0.0
            𝑀[i,j], Θ[i,j] = _rateofchange(tmp)
        end
    end
end

change = SimpleSDMResponse(𝑀, A)
replace!(change, 0.0 => nothing)

qc = rescale(change, collect(0.0:0.01:1.0))

plot(change, frame=:box, title="Rate of change", clim=Tuple(quantile(collect(change), [0.1,0.90])), c=:roma, dpi=600)

plot(qc, c=:roma, dpi=600, title="Quantiles of the rate of change")