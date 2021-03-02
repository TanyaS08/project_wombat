

#=

    We expect two types of layouts of the data
    which determines the wombling method

=#

abstract type GridStructure end; #New Parent type defining relationship between co-ords

abstract type Regular <: GridStructure end; #'Perfect' grid
abstract type Random <: GridStructure end; #Points 'random'

const AG = ArchGDAL

#1 test layout of co-ords to see how they are arranged

#2 Assign and return result as type of GridStructure i.e. Lattice ot Triangular



#

## STEP 1: Create the network surface

## STEP 2: Connecting the networks in space

# TODO how to Delaunay Triangulate

#sort sites into buckets
#sort along 1 axis and if equal then by the Second
sort(amphdata.Lat, by=sort(amphdata.Long))

amph_coord = convert(Matrix, amphdata[1:10, [:Lat, :Long]]);

mesh = delaunay(amph_coord);

mesh.neighbors
mesh.points
mesh.simplices

C = [amph_coord[i,1] for i in mesh.simplices];
D = [amph_coord[i,2] for i in mesh.simplices];
#amphdata[1:10, :Salamandra_salamandra]
#Z = [amph_coord[i,2] for i in mesh.simplices];

mesh.vertex_to_simplex

using Makie
color = rand(size(mesh.points, 1))
scene = Makie.mesh(mesh.points, mesh.simplices, color=color, shading=false, scale_plot=false)
Makie.wireframe!(scene[end][1], color=(:black, 0.6), linewidth=5)

#triangualte cells

#merge cells into rows



# For lattice - still need to evaluate which points are 'realted' to which

## STEP 3: Traverse and 'collate' the surface

amphdata.Long

#Range Scale between 0 & 1
x_scaled = (amphdata.Long .- minimum(amphdata.Long))/(maximum(amphdata.Long) - minimum(amphdata.Long));
y_scaled = (amphdata.Lat .- minimum(amphdata.Lat))/(maximum(amphdata.Lat) - minimum(amphdata.Lat));

unique(sort(x_scaled))
unique(y_scaled)

minimum(x_scaled)

struct Rate
    𝑚::Float32
    θ::Float32
end

A = [1 2 3 4; 5 6 7 8; 9 10 11 12]

Assemblage(occ = amphdata)
SpatialEcology.parsesingleDataFrame(amphdata)

filter(:Salamandra_salamandra => x -> x > 0, amphdata)

histogram2d(amphdata.Long[1:10, :],
amphdata.Lat[1:10, :], c=:viridis)


## WOMBLE!

"""
    TODO - This calulates the rate of change (𝑚) for 3 points
"""

#C::A is the co-ordiantes
#Z::X is the Z values
#TODO Add type to differentiate between co-ord layout i.e. facilitate multiple dispatch

C = cat(C[1,:],D[1,:],[1,1,1], dims =(2, 2));
Z = rand((0:10), (10,3));

RateOfChange.(C,D,Z)

function RateOfChange(Lat::Vector, Long::Vector, Z::Vector) #where {Random <: GridStructure, X <: Vector{Float64}}

    C = cat(Lat,Long,[1,1,1], dims =(2, 2))

    coeff = Base.inv(C) * Z;
    𝑋 = sum(C[:,1])/3; #X co-ord
    𝑌 = sum(C[:,2])/3; #Y co-ord

    ∂𝑋 = coeff[2]*𝑌 + coeff[3]
    ∂𝑌 = coeff[1]*𝑋 + coeff[3]

    𝑚 = sqrt(exp2(∂𝑋) + exp2(∂𝑌))

    if ∂𝑋 < 0
        Δ = 180;
    else
        Δ = 0;
    end
    θ = atan(∂𝑋/∂𝑌) + Δ;
    return [𝑚, θ]
end

"""
    TODO - This calulates the rate of change (𝑚) for a lattice
"""
function RateOfChange(A::Regular) where {Regular <: GridStructure}

#in theory this traverses the entire matrix
    for i in 1:(size(A, 2) - 1) #col
        for j in 1:(size(A, 1) - 1) #row

            𝑍₄ = A[j, i];
            𝑍₃ = A[j, i + 1];
            𝑍₂ = A[j + 1, i + 1];
            𝑍₁ = A[j + 1, i];

            ∂𝑋 = 𝑍₂ - 𝑍₁ + 0.5(𝑍₁ - 𝑍₂ + 𝑍₃ - 𝑍₄);
            ∂𝑌 = 𝑍₄ - 𝑍₁ + 0.5(𝑍₁ - 𝑍₂ + 𝑍₃ - 𝑍₄);

            𝑚 = sqrt(exp2(∂𝑋) + exp2(∂𝑌));

            if ∂𝑋 < 0
                Δ = 180;
            else
                Δ = 0;
            end
            return Δ

            θ = atan(∂𝑋/∂𝑌) + Δ,

            push!(Rate(𝑚, θ))

        end
    end
end


## STEP 5: Threshold values

## STEP 6: Signif of candidiate boundaries

## STEP 7: Test if boundaries are connected
