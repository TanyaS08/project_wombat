

#=

    We expect two types of layouts of the data
    which determines the wombling method

=#

const AG = ArchGDAL

#1 test layout of co-ords to see how they are arranged

#2 Assign and return result as type of GridStructure i.e. Lattice ot Triangular



#

## STEP 1: Create the network surface

## STEP 2: Connecting the networks in space

# TODO how to Delaunay Triangulate - port from SciPy??

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

#Range Scale between 0 & 1
x_scaled = (amphdata.Long .- minimum(amphdata.Long))/(maximum(amphdata.Long) - minimum(amphdata.Long));
y_scaled = (amphdata.Lat .- minimum(amphdata.Lat))/(maximum(amphdata.Lat) - minimum(amphdata.Lat));

unique(sort(x_scaled))
unique(y_scaled)

minimum(x_scaled)

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

[RateOfChange(C[i,:], D[i,:], Z[i,:]) for i in 1:10]

## STEP 5: Threshold values

## STEP 6: Signif of candidiate boundaries

## STEP 7: Test if boundaries are connected
