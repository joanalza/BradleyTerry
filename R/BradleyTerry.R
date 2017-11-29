#' An S4 class to represent distributions based on Bradley Terry distribution model.
#'
#' @slot abilities The logarithm of weights parameter of the learned population.
#'

setClass(
  Class = "BradleyTerry",
  representation = representation(abilities = "numeric", sampling = "character", indices = "data.frame", burnIn = "numeric")
)
# GENERIC METHODS ---------------------------------------------------------
setMethod(
  f = "simulate",
  signature = "BradleyTerry",
  definition = function(object, nsim = 1, seed = NULL, ...) {
    switch(
      object@sampling,
      Random = {
        abilities <- object@abilities # create a copy to manipulate it
        new.population <- lapply(
          1:nsim, FUN = function(x) {
            randomPermutation(length(abilities))
          }
        )
      },
      Metropolis = {
        abilities <- object@abilities
        n <- length(abilities)
        new.population <- list()
        for( i in 1:nsim){
          sigma0 <- randomPermutation(n)
          aux0 <- apply(object@indices, MARGIN = 1, FUN = calcProb, sigma = sigma0, ability = abilities)
          logprob0 <- sum(aux0)
          for ( j in 1:object@burnIn){
            #browser()
            ind <- sample(1:n, 2, replace = FALSE)
            sigma1 <- swap(sigma0, ind[1], ind[2])
            aux1 <- apply(object@indices, MARGIN = 1, FUN = calcProb, sigma = sigma1, ability = abilities)
            logprob1 <- sum(aux1)
            probRatio <- exp(logprob1 - logprob0)
            if(runif(1) < probRatio){
              sigma0 <- sigma1
              aux0 <- apply(object@indices, MARGIN = 1, FUN = calcProb, sigma = sigma0, ability = abilities)
              logprob0 <- sum(aux0)
            }
          }
          new.population <- append(new.population, sigma0)
        }
      },
      Heuristic = {
        #browser()
        abilities <- object@abilities
        n <- length(abilities)
        new.population <- list()
        for (i in 1:nsim){
          order.vec <- sample(1:n)                # Get random order of n elements.
          sigma <- c(head(order.vec,1))
          #browser()
          for (j in 2:n){
            bool.vec <- c()
            enter.bool <- FALSE
            new.elem <- order.vec[j]
            index <- sample(length(sigma),1)
            while( !enter.bool ){      # Loop infinituak ekiditeko.
              prob <- exp(abilities[new.elem]) / ( exp(abilities[new.elem]) + exp(abilities[sigma[index]]) )
              if ( runif(1) >= prob){            # Atzetik jarri behar
                bool.vec[index] <- 2    #Atzetik aldagaia
                index <- index + 1
              }else{
                bool.vec[index] <- 1    #Aurretik aldagaia
                index <- index - 1
              }
              if (index == 0){
                sigma <- c(new.elem,sigma)
                enter.bool <- TRUE
              }else if (index > length(sigma)){
                sigma <- c(sigma,new.elem)
                enter.bool <- TRUE
              }else if ( !is.na(bool.vec[index]) ){
                if(bool.vec[index] == 1){
                  sigma <- c(sigma[1:index - 1],c(new.elem,sigma[index:length(sigma)]))
                }else{
                  index <- index + 1
                  sigma <- c(sigma[1:index - 1],c(new.elem,sigma[index:length(sigma)]))
                }
                enter.bool <- TRUE
              }
            }
          }
          new.population <- append(new.population, permutation(sigma))
        }
      }
    )
    return(new.population)
  }
)

# CONSTRUCTOR -------------------------------------------------------------

#' This function creates an object of class \code{\linkS4class{BradleyTerry}}
#'
#' @family EDA
#' @param data The dataframe containing the initial-population to contruct the Bradley Terry model.
#' @param ... Ignored
#' @return An object of class \code{\linkS4class{BradleyTerry}} that includes the weight parameter vector of the given data
#'
bradleyTerry <- function(data, burnIn = 100, sampling = "Heuristic", ...) {
  if (class(data) == "list") {
    P <- length(data[[1]]@permutation) # Problem size
    N <- length(data) # Pop size
    X <- P - 1 # Pop size
    Y <- 2 #Problem size
    
    player1 <- c()
    player2 <- c()
    win1 <- c()
    win2 <- c()
    
    for (i in 1:X) {
      for (j in Y:P) {
        player1 <- c(player1,i)
        player2 <- c(player2,j)
        val <- checkProb(data, i,j)
        win1 <- c(win1, val)
        win2 <- c(win2, N - val)
      }
      Y <- Y + 1
    }
    d <- data.frame(player1,player2,win1,win2)
    d$player1 <- factor(d$player1, levels = 1:P)
    d$player2 <- factor(d$player2, levels = 1:P)
    model <-
      BTm(
        cbind(win1,win2), player1 = player1, player2 = player2, formula = ~ player, id =
          "player", data = d
      )
    obj <-
      new("BradleyTerry", abilities = BTabilities(model)[,1], sampling = sampling, indices = d[,1:2], burnIn = burnIn)
  }else{
    stop("The data must be a list")
  }
  return (obj)
}

#' How many times object 1 is earlier than object 2 in a permutation
#' @param x numeric list (permutation).
#' @param obj1 number of the list.
#' @param obj2 number of the list.
checkProb <- function(x, obj1, obj2) {
  sum <- 0
  N <- length(x)
  for (i in 1:N) {
    elem1 <- which(as.numeric(x[[i]]) == obj1)
    elem2 <- which(as.numeric(x[[i]]) == obj2)
    if (elem1 < elem2) {
      sum <- sum + 1
    }
  }
  return(sum)
}

#' Calculates probabily of two elements
#' @param vector with two indices
#' @param permutation
#' @param probabilities of the elements
calcProb <- function (ind, sigma, ability){
  i <- as.numeric(ind[1])
  j <- as.numeric(ind[2])
  w_sigma_i <- exp(ability[sigma[i]])
  w_sigma_j <- exp(ability[sigma[j]])
  return ( log(w_sigma_i / (w_sigma_i + w_sigma_j)) )
}

#' Calculates probability of a permutation
#' @param model to take abilities
#' @param permutation
#' 
#' @return probability
getProbability <- function (model, sigma){
  aux0 <- apply(model@indices, MARGIN = 1, FUN = calcProb, sigma = sigma, ability = model@abilities)
  return (exp(sum(aux0)))
}