## Function to create a set of grid cells based on a set of x and y coordinates
#
#	x				-	Vector of x-coordinates of the grid cells
#	y				-	Vector of y-coordinates of the grid cells
#
MakeGridFromCoords <- function(x, y)
{
  out <- .External(
    "MakeGridFromCoords",
    x = as.double(x),
    y = as.double(y),
    PACKAGE = "ecomodtools")
  data.frame(
    x1 = out$x1,
    x2 = out$x2,
    y1 = out$y1,
    y2 = out$y2)
}

# Create a data frame of cell boundary coordinates of cell of a 5x5
# lattice with regular cell length of one unit
x.coords <- seq(0.0, 5.0, 1.0)
y.coords <- seq(0.0, 5.0, 1.0)
MakeGridFromCoords(x.coords, y.coords)

# #A data frame with (length(x.coords) - 1) * (length(y.coords) - 1) rows containing the boundary coordinates for each cell in the lattice:
#   
#   x1	
# Lower x-coordinate of each cell in the lattice
# x2	
# Upper x-coordinate of each cell in the lattice
# y1	
# Lower y-coordinate of each cell in the lattice
# y2	
# Upper y-coordinate of each cell in the lattice
#   