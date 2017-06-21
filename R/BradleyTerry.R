#' An S4 class to represent distributions based on Bradley Terry distribution model.
#'
#' @slot parameters The weight parameters of the learned population.
#' 

setClass(
  Class="BradleyTerry", 
  representation=representation(parameters="numeric")
)
# GENERIC METHODS ---------------------------------------------------------
setMethod(
  f="simulate", 
  signature="BradleyTerry", 
  definition=function(object, nsim=1, seed=NULL, ...) { # Random Permutation
    new.population <- c()
    parameters <- object@parameters # create a copy to manipulate it
    new.population <- randomPermutation(length(parameters))
    return(new.population)
  })

# CONSTRUCTOR -------------------------------------------------------------

#' This function creates an object of class \code{\linkS4class{BradleyTerry}}
#' 
#' @family EDA
#' @param data The dataframe containing the initial-population to contruct the Bradley Terry model.
#' @param ... Ignored
#' @return An object of class \code{\linkS4class{BradleyTerry}} that includes the weight parameter vector of the given data
#' 
bradleyTerry <- function(data, maxIter, ...) {
  if(missing(maxIter)){
    maxIter <- 10000
  }
  if(class(data) == "list"){
    #data <- c(data, permutation(rev(data[[1]]@permutation)))
    # We can create the object
    # Transform data to a matrix
    P <- length(data[[1]]@permutation) # Problem size
    N <- length(data) # Pop size
    X <- P-1 # Pop size
    Y <- 2 #Problem size
    
    player1 <- c()
    player2 <- c()
    win1 <- c()
    win2 <- c()
    
    for(i in 1:X){
      for(j in Y:P){
        player1 <- c(player1,i)
        player2 <- c(player2,j)
        val <- checkProb(data, i,j)
        win1 <- c(win1, val)
        win2 <- c(win2, N - val)
      }
      Y <- Y + 1
    }
    return(data.frame(player1,player2,win1,win2))
  }else{
    stop("The data must be a list")
  }
}

checkProb <- function(x, obj1, obj2) {
  sum <- 0
  N <- length(x)
  for (i in 1:N){
    elem1 <- which(as.numeric(x[[i]]) == obj1)
    elem2 <- which(as.numeric(x[[i]]) == obj2)
    if (elem1 < elem2){
      sum <- sum +1
    }
  }
  return(sum)
}