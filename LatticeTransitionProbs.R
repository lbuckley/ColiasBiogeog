# SEE: https://rdrr.io/rforge/ecomodtools/man/LatticeTransitionProbs.html

LatticeTransitionProbs(x1, x2, y1, y2, func, approx.method = "AA", boundary = "absorbing", max.prob = .Machine$double.eps^0.25, initial.step = 10000, step.size = 10000, max.rel.var = .Machine$double.eps^0.25, max.runs = 1000000, params = NULL, ...)

## Function to calculate cell-to-cell transition probabilities
#
#	x1				-	Lower x-coordinates of lattice cells
#	x2				-	Upper x-coordinates of lattice cells
#	y1				-	Lower y-coordinates of lattice cells
#	y2				-	Upper y-coordinates of lattice cells
#	func				-	Dispersal function can be a character vector with the following values
#							"gaussian" : Gaussian dispersal function
#	approx.method			-	Approximation method to use
#							"CC" : Centroid-to-centroid dispersal
#							"CA" : Centroid-to-area dispersal
#							"AC" : Area-to-centroid dispersal
#							"AA" : Area-to-area dispersal
#	max.prob			-	Stopping condition for the calculation of infinite series
#	initial.step			-	MC integration parameter - number of runs to do in the first step
#	step.size			-	MC integration parameter - number of runs to do in one time step
#	max.rel.var			-	MC integration parameter - maximum relative variance before stopping MC integration
#	max.runs			-	MC integration parameter - maximum number of runs before stopping MC integration
#	params				-	Parameter vector used in dispersal calculation
#	...				-	Extra parameters to pass to the dispersal function
#
LatticeTransitionProbs <- function(x1, x2, y1, y2, func, approx.method = "AA", boundary = "absorbing", max.prob = .Machine$double.eps^0.25,
                                   initial.step = 10000, step.size = 10000, max.rel.var = .Machine$double.eps^0.25, max.runs = 1000000, params = NULL, ...)
{
  out <- NULL
  # Check to see if function input is a character vector
  if(is.character(func))
  {
    if(tolower(func[1]) == "gaussian")
    {
      # Function is Gaussian - call relevant C++ function
      out <- t(.External("GaussianLatticeTransition",
                         x1 = as.double(x1),
                         x2 = as.double(x2),
                         y1 = as.double(y1),
                         y2 = as.double(y2),
                         alpha = as.double(params)[1],
                         approx.method = toupper(as.character(approx.method)[1]),
                         boundary = tolower(as.character(boundary)[1]),
                         max.prob = as.double(max.prob)[1],
                         PACKAGE = "ecomodtools"))
    }
    else
    {
      stop("unsupported function name")
    }
  }
  # ... otherwise resort to Monte Carlo integration
  else
  {
    # Function wrapper for centroid-to-centroid dispersal
    CCDisp <- function(from.x1, from.x2, from.y1, from.y2,
                       to.x1, to.x2, to.y1, to.y2,
                       func.in, params, initial.step, step.size, max.rel.var, max.runs, ...)
    {
      func.in(
        from = c((from.x1 + from.x2) / 2.0, (from.y1 + from.y2) / 2.0),
        to = c((to.x1 + to.x2) / 2.0, (to.y1 + to.y2) / 2.0),
        params = params, ...)
    }
    # Function wrapper for centroid-to-area dispersal
    CADisp <- function(from.x1, from.x2, from.y1, from.y2,
                       to.x1, to.x2, to.y1, to.y2,
                       func.in, params, initial.step, step.size, max.rel.var, max.runs, ...)
    {
      # Wrapper function for centroid-to-area integration
      CAWrap <- function(x, params, from, ...)
      {
        func.in(
          from = from,
          to = c(x[1], x[2]),
          params = params, ...)
      }
      # Perform the Monte Carlo integration
      MCIntegration(func = CAWrap, lower = c(to.x1, to.y1), upper = c(to.x2, to.y2),
                    initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs,
                    params = params, from = c((from.x1 + from.x2) / 2.0, (from.y1 + from.y2) / 2.0), ...)$val
      
    }
    # Function wrapper for area-to-centroid dispersal
    ACDisp <- function(from.x1, from.x2, from.y1, from.y2,
                       to.x1, to.x2, to.y1, to.y2,
                       func.in, params, initial.step, step.size, max.rel.var, max.runs, ...)
    {
      # Wrapper function for area-to-centroid integration
      ACWrap <- function(x, params, to, ...)
      {
        func.in(
          from = c(x[1], x[2]),
          to = to,
          params = params, ...)
      }
      # Perform the Monte Carlo integration
      MCIntegration(func = ACWrap, lower = c(from.x1, from.y1), upper = c(from.x2, from.y2),
                    initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs,
                    params = params, to = c((to.x1 + to.x2) / 2.0, (to.y1 + to.y2) / 2.0), ...)$val
    }
    # Function wrapper for area-to-area dispersal
    AADisp <- function(from.x1, from.x2, from.y1, from.y2,
                       to.x1, to.x2, to.y1, to.y2,
                       func.in, params, initial.step, step.size, max.rel.var, max.runs, ...)
    {
      # Wrapper function for area-to-area integration
      AAWrap <- function(x, params, ...)
      {
        func.in(
          from = c(x[1], x[2]),
          to = c(x[3], x[4]),
          params = params, ...)
      }
      # Perform the Monte Carlo integration
      MCIntegration(func = AAWrap, lower = c(from.x1, from.y1, to.x1, to.y1), upper = c(from.x2, from.y2, to.x2, to.y2),
                    initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs,
                    params = params, ...)$val
    }
    # List of dispersal functions
    disp.funcs <- list(CC = CCDisp, CA = CADisp, AC = ACDisp, AA = AADisp)
    # Calculate the cell-to-cell transition probabilities using Monte Carlo simulation
    out <- t(.External("LatticeTransition",
                       x1 = as.double(x1),
                       x2 = as.double(x2),
                       y1 = as.double(y1),
                       y2 = as.double(y2),
                       func = as.function(func),
                       approx.method = toupper(as.character(approx.method)[1]),
                       boundary = tolower(as.character(boundary)[1]),
                       max.prob = as.double(max.prob)[1],
                       initial.step = as.integer(initial.step)[1],
                       step.size = as.integer(step.size)[1],
                       max.rel.var = as.double(max.rel.var)[1],
                       max.runs = as.integer(max.runs)[1],
                       params = params,
                       disp.funcs = disp.funcs,
                       PACKAGE = "ecomodtools", ...))
  }
  out
}