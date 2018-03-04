#' An S4 class to represent distributions based on Bradley Terry distribution model.
#'
#' @slot abilities The logarithm of weights parameter of the learned population.
#'

setClass(
  Class = "MixedPLBT",
  representation = representation(
    theta = "numeric", PL = "PlackettLuce", BT = "BradleyTerry"
  )
)
# GENERIC METHODS ---------------------------------------------------------
setMethod(
  f = "simulate",
  signature = "MixedPLBT",
  definition = function(object, nsim = 1, seed = NULL, ...) {
    nBin.BT <- rbinom(1, nsim, object@theta)
    nBin.PL <- nsim - nBin.BT
    if (nBin.BT == 0){
      BT.simulate <- NULL
    }else{
      BT.simulate <- simulate(object@BT, nBin.BT)
    }
    if (nBin.PL == 0){
      PL.simulate <- NULL
    }else{
      PL.simulate <- simulate(object@PL, nBin.PL)
    }
    new.population <- append( BT.simulate, PL.simulate)
    return (new.population)
  }
)


# CONSTRUCTOR -------------------------------------------------------------

#' This function creates an object of class \code{\linkS4class{BradleyTerry}}
#'
#' @family EDA
#' @param data The dataframe containing the initial-population to contruct the Bradley Terry and Plackett-Luce models.
#' @param ... Ignored
#'
mixedPLBT <- function(data, file, ...) {
  if (class(data) == "list") {
    d <- length(data[[1]]@permutation) # Problem size
    M <- length(data) # Pop size
    Vr <- ((d^2) -1)/12
    
    data.matrix <- t(sapply(data, '[', 1:max(sapply(data, length)))) # Population list as matrix
    Xf.var <- apply(data.matrix, MARGIN = 2, var)
    theta <- 1 - (mean(Xf.var)/ Vr)    # Nogueira's index
    theta <- max(0, theta)
    
    f <- file(description = file, open = "a")
    write.table(theta, file = f, col.names = FALSE, sep = ",")
    close(f)
    
    BT <- bradleyTerry(unique(data))
    PL <- plackettLuce(data, maxIter = 1000)
    obj <-
      new("MixedPLBT", theta = theta, BT = BT, PL = PL)
  }else{
    stop("The data must be a list")
  }
  return (obj)
}