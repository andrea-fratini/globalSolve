# Define the SmolyakKernel class
SmolyakKernel <- setRefClass(
  "SmolyakKernel",
  fields = list(
    D = "integer",          # Dimensions
    mu = "integer",         # Index of mu
    N = "integer",          # Number of Grid Points
    xbnds = "list",         # Bounds of dimensions of x
    zbnds = "list",         # Bounds of dimensions of z
    GridIdx = "list",       # Input for grid construction
    BasisIdx = "list"       # Input to construct Basis Funs for set of grid points
  )
)

# Constructor functions
SmolyakKernel$methods(
  initialize = function(D, mu_level, xbnds = NULL, zbnds = NULL, HD = FALSE) {
    if (missing(D) || missing(mu_level)) {
      stop("D and mu_level must be provided")
    }
    
    muvec <- rep(mu_level, D)
    
    if (HD) {
      GridIdx <- HDSmolIdx(muvec)
      NGP <- NumGridPts(GridIdx, muvec)
    } else {
      result <- SmolIdx(muvec)
      NGP <- result$NumGridPts
      GridIdx <- result$iList
    }
    
    if (is.null(xbnds)) {
      xbnds <- rep(list(c(-1.0, 1.0)), D)
    }
    
    if (is.null(zbnds)) {
      zbnds <- rep(list(c(-1.0, 1.0)), D)
    }
    
    BasisIdx <- lapply(1:NGP, function(r) integer(D))
    BasisIdx <- BasisIdxUpdate(BasisIdx, GridIdx, muvec)
    
    .self$D <- D
    .self$mu <- muvec
    .self$N <- NGP
    .self$xbnds <- xbnds
    .self$zbnds <- zbnds
    .self$GridIdx <- GridIdx
    .self$BasisIdx <- BasisIdx
  }
)

# Helper functions
chebtheta <- function(n) {
  if (n == 1) return(0.5 * pi)
  return((0:(n-1)) * pi / (n-1))
}

chebnodes <- function(n) {
  return(round(cos(chebtheta(n)), digits = 14))
}

m_i <- function(i) {
  if (i == 1) return(1)
  return(2^(i-1) + 1)
}

newpts <- function(i) {
  if (i <= 2) return(i)
  return(2^(i-2))
}

grid_A_i <- function(i) {
  if (i == 1) return(0.0)
  if (i == 2) return(chebnodes(m_i(i))[seq(1, m_i(i), by = 2)])
  return(chebnodes(m_i(i))[seq(2, m_i(i), by = 2)])
}

# Calculate Number of Grid Points & Construct Indices of Smolyak Grid
SmolIdx <- function(mu) {
  max_mu <- max(mu)
  iList <- list()
  ibar <- integer()
  n <- integer()
  NumGridPts <- 0
  sum_i <- 0
  
  for (i in 1:length(mu)) {
    sum_i <- sum_i + i
    if (sum_i > length(mu) + max_mu) {
      ibar <- integer()
      n <- integer()
      sum_i <- 0
      break
    }
    ibar <- c(ibar, i)
    iList <- append(iList, list(ibar))
    NumGridPts <- NumGridPts + prod(n)
    ibar <- integer()
    n <- integer()
    sum_i <- 0
  }
  
  return(list(NumGridPts = NumGridPts, iList = iList))
}

# Disjoint sets that define indexes for creation of Basis Indices
A_pidx <- function(ibar) {
  A <- list()
  lb <- c(1, 2)
  lb <- c(lb, 2^(3:ibar - 1) - 2^(3:ibar - 2) + 2)
  ub <- c(1)
  ub <- c(ub, 2^(2:ibar - 1) + 1)
  for (j in 1:ibar) {
    A <- append(A, list(seq(lb[j], ub[j])))
  }
  return(A)
}

# Make Basis Function Indexes
BasisIdxUpdate <- function(Binds, GridIdx, mu) {
  max_mu <- max(mu)
  A <- A_pidx(max_mu + 1)
  zT <- 0
  
  for (i in 1:length(GridIdx)) {
    for (j in expand.grid(lapply(GridIdx[[i]], function(k) A[[k]]))) {
      zT <- zT + 1
      for (w in 1:length(mu)) {
        Binds[[zT]][w] <- j[[w]] - 1
      }
    }
  }
  
  return(Binds)
}

SmolyakKernel()

# Define the SmolyakGrid class
SmolyakGrid <- setRefClass(
  "SmolyakGrid",
  fields = list(
    sk = "SmolyakKernel",  # SmolyakKernel
    N = "integer",         # Number of Grid Points
    grid = "list"          # Smolyak Grid on z
  ),
  methods = list(
    initialize = function(sk = NULL, D = 2, mu_level = 3, mu = NULL, xbnds = NULL, HD = FALSE) {
      if (!is.null(sk)) {
        # Call with Smolyak Kernel setup
        .self$sk <- sk
        .self$N <- NumGridPts(sk$GridIdx, sk$mu)
        .self$grid <- lapply(1:.self$N, function(r) numeric(sk$D))
        gridUpdate(.self$grid, sk$GridIdx, sk$mu)
      } else if (!is.null(D) && !is.null(mu_level)) {
        # Set up Smolyak Kernel internally
        .self$sk <- SmolyakKernel(D = D, mu_level = mu_level, HD = HD)
        .self$N <- NumGridPts(.self$sk$GridIdx, .self$sk$mu)
        .self$grid <- lapply(1:.self$N, function(r) numeric(.self$sk$D))
        gridUpdate(.self$grid, .self$sk$GridIdx, .self$sk$mu)
      } else if (!is.null(mu)) {
        # Set up Smolyak Kernel internally
        .self$sk <- SmolyakKernel(mu = mu, HD = HD)
        .self$N <- NumGridPts(.self$sk$GridIdx, .self$sk$mu)
        .self$grid <- lapply(1:.self$N, function(r) numeric(.self$sk$D))
        gridUpdate(.self$grid, .self$sk$GridIdx, .self$sk$mu)
      } else if (!is.null(D) && !is.null(mu_level) && !is.null(xbnds)) {
        # Set up Smolyak Kernel internally with xbnds
        .self$sk <- SmolyakKernel(D = D, mu_level = mu_level, xbnds = xbnds, HD = HD)
        .self$N <- NumGridPts(.self$sk$GridIdx, .self$sk$mu)
        .self$grid <- lapply(1:.self$N, function(r) numeric(.self$sk$D))
        gridUpdate(.self$grid, .self$sk$GridIdx, .self$sk$mu)
      } else if (!is.null(mu) && !is.null(xbnds)) {
        # Set up Smolyak Kernel internally with xbnds
        .self$sk <- SmolyakKernel(mu = mu, xbnds = xbnds, HD = HD)
        .self$N <- NumGridPts(.self$sk$GridIdx, .self$sk$mu)
        .self$grid <- lapply(1:.self$N, function(r) numeric(.self$sk$D))
        gridUpdate(.self$grid, .self$sk$GridIdx, .self$sk$mu)
      } else {
        stop("Invalid arguments provided. Either 'sk' or 'D' and 'mu_level' must be provided.")
      }
    }
  )
)

# Helper function to update the grid
gridUpdate <- function(H, inds, mu) {
  zH <- 0
  for (i in 1:length(inds)) {
    for (j in expand.grid(lapply(inds[[i]], function(k) grid_A_i(k)))) {
      zH <- zH + 1
      for (w in 1:length(mu))) {
        H[[zH]][w] <- j[[w]]
      }
    }
  }
return(H)
}

# Helper function to calculate the number of grid points
NumGridPts <- function(GridIdx, mu) {
  return(length(GridIdx))
}

# Example usage
# Create a SmolyakGrid object with D and mu_level
sg <- SmolyakGrid(D = 2, mu_level = 3)
print(sg$N)  # Number of grid points
print(sg$grid)  # Smolyak grid