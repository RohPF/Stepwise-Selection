library(lmtest)

#' This function gives the user easy access to p-values for specific predictors from a linear model.  It is mainly used
#' to easily pass p-values to other functions.
#'
#' @param pred the predictor from the linear model whose p-value is desired
#' @param fit a linear model of type "glm" containing the desired predictor
#' @param step require choosing ether step forward or backward
#' @return the p-value of our interest according to the step

extractp <- function(pred, fitCurrent, step = 'fwd') {

  # Foward with Rao Score
  if (step == 'fwd'){
    before = fitCurrent
    after = fMaker(pred, fitCurrent, add=T) # Add
    pvalue = anova(before, after, test='Rao')[2,6] # Rao
    names(pvalue) = pred
    return(as.data.frame(pvalue))
  }
  
  # Backward with Wald Test
  else if (step == 'bwd')
    before = fitCurrent
    after = fMaker(pred, fitCurrent, add=F) # Remove
    pvalue <- lmtest::waldtest(before, after)[2,4] # Wald
    names(pvalue) = pred
    return(pvalue)
}

#' Create formulas from strings
#' @param pred the predictor to be added or removed from the current model
#' @param fitCurrent the current model to be updated
#' @param add by default adds the predictor to to the model.  add=F removes the predictor from the model
#' @return the updated model of type "glm"

fMaker <- function(pred, fitCurrent, add=T) {
  addNew <- as.formula(paste(".~.+", pred))
  remNew <- as.formula(paste(".~.-", pred))
  if(add) {
    return(update(fitCurrent, addNew))
  }  else {
    return(update(fitCurrent, remNew))
  }
}

#' Adds a single predictor to a linear model based on its p-value
#' The comparison is done using rao score test
#' This function will try and add a new predictor to a current model.  A predictor will be added if it has minimum p-value among all predictors and its p-value is below a certain threshold
#'
#' @param fitCurrent the current model of type "glm"
#' @param fullmodel a linear model containing all possible predictors.  Typically of the form glm(y~., data=data)
#' @param aEnter the threshold for adding the predictor, set to 0.1 be default
#' @param forcedOut vector of predictors that will be forced out of the final model
#' @return an updated linear model of type "glm"

stepfwd <- function(fitCurrent, fullmodel, aEnter = 0.1, forcedOut = NULL) {
  predsInModel <- rownames(anova(fitCurrent))                                                   #list of predictors in current model
  predsFull <- rownames(anova(fullmodel))                                                       #list of predictors in full model
  predsNotInModel <- setdiff(predsFull, predsInModel)                                           #list of predictors not in current model, not including forced out predictors
  pvalues <- unlist(sapply(predsNotInModel, function(x) as.numeric(extractp(x, fitCurrent, step = 'fwd'))))   #takes each predictor not in the current model, creates a new lm which includes it, and stores its respective p-value. Really ugly, but the only way I could make it work.
  print(as.data.frame(round(pvalues,5)))
  cat('\n')
  if(length(pvalues)==0) return(fitCurrent)
  toAdd <- pvalues[which(pvalues==min(pvalues, na.rm = TRUE))]                                         #possible new predictor
  if(as.numeric(toAdd) <= aEnter) {
    cat("+++++ Add predictor", names(toAdd), "+++++", "\n")
    print(summary(fMaker(names(toAdd), fitCurrent))$coefficients, digits = 4)
    cat("\n")
    return(fMaker(names(toAdd), fitCurrent))                                                    #updates and returns new model with additional predictor
  }
  return(fitCurrent)
}

#' Removes a single predictor from a linear model based on its p-value
#'
#' This function will try and remove a single predictor from a current linear model.  A predictor will be removed if it has maximal p-value and its p-value is greater than a certain threshold.
#'
#' @param fitCurrent the current model of type "lm"
#' @param fullmodel a linear model containing all possible predictors.  Typically of the form lm(y~., data=data)
#' @param forcedIn vector of predictors that will be forced into the final model
#' @param aRemove the threshold for removing the predictor, set to 0.1 by default
#' @return an updated linear model of type "lm"

stepbwd <- function(fitCurrent, fullmodel, aRemove = 0.1, forcedIn = NULL) {
  predsIncluded <- rownames(anova(fitCurrent))                                               #predictors in current model
  predsIncluded <- predsIncluded[(predsIncluded != "NULL")]                                  #removes "Null" from predictors
  predsIncluded <- setdiff(predsIncluded, intersect(predsIncluded, forcedIn))                #makes sure no forced in predictors get removed
  pvalues <- unlist(sapply(predsIncluded, function(x) as.numeric(extractp(x, fitCurrent, step = 'bwd'))))  #checks the p-value for each predictor in current model
  print(as.data.frame(round(pvalues,5)))
  cat('\n')
  if(length(pvalues)==0) return(fitCurrent)                                                  #returns current model if there are no more possible predictors to add
  toRemove <- pvalues[which(pvalues==max(pvalues, na.rm = TRUE))]                                     #selects the predictor with maximal p-value
  if(length(toRemove)==0) return(fitCurrent)
  if(as.numeric(toRemove) > aRemove){
    cat("----- Remove predictor", names(toRemove), "-----", "\n")
    print(summary(fMaker(names(toRemove), fitCurrent, add=FALSE))$coefficients, digits = 4)
    cat("\n")
    return(fMaker(names(toRemove), fitCurrent, add=FALSE))                                   #returns an updated model if the p-value is above the threshold
  }
  return(fitCurrent)                                                                       #else, returns original model
}

#' Selects the best predictors for a linear model based on p-values
#'
#' This function will attempt to create the "best' linear model by finding the most significant predictors. A predictor will be included/excluded in the final model if when it is added/removed its p-value is below/above a certain threshold.
#' @usage pStepwise(response, fullmodel, aEnter = 0.1,
#'                  aRemove = 0.1, forcedIn = NULL, forcedOut = NULL,
#'                  method = "both", plotRes = FALSE)
#' @param response the response variable of interest in the model
#' @param fullmodel a linear model containing all possible predictors, typically of the form lm(y ~ ., data = data)
#' @param aEnter the threshold for adding new predictors, set to 0.1 by default
#' @param aRemove the threshold for removing predictors from the current model, set to 0.1 by default
#' @param forcedIn a vector of predictors that will be forced into the final model regardless of their p-values
#' @param forcedOut a vector of predictors that will not be included in the final model regardless of their p-values
#' @param method "forward" will only add predictors, "backward" will only remove predictors and "both" will perfrom stepwise.  "both" by default
#' @param plotRes option to print residual plots.  plotRes = TRUE will include normailty plot and fit vs. residual plots. FALSE by defult.
#' @import stats
#' @return a linear model of type lm containing the "best" predictors
#'
pStepwise <- function(response, initmodel, fullmodel, aEnter = 0.3, aRemove = 0.35,
                      forcedIn = NULL, forcedOut = NULL, method = "stepwise") {
  
  fitBwd <- initmodel                           #creates an empty model to begin with (poor name choice but it makes the while loop easier)
  for(pred in forcedIn) fitBwd <- fMaker(pred, fitBwd)                             #adds in forced predictors to initial model
  for(pred in forcedOut) fullmodel <- fMaker(pred, fullmodel, add=F)               #removes forced out predictors from fullmodel
  if(method=="forward") aRemove = 1                                                #makes it impossible to remove predictors
  if(method=="backward") {                                                         #section for backward selection
    fitFwd <- fullmodel
    while(TRUE) {
      fitBwd <- stepbwd(fitFwd, fullmodel, aRemove = aRemove, forcedIn = forcedIn)
      if(identical(fitFwd, fitBwd)==T) {
        cat("======== Final model ========", "\n")
        print(fitFwd$call)
        cat("Predictors forced in: ", forcedIn, "\n")
        cat("Predictors forced out: ", forcedOut, "\n", "\n")
        print(summary(fitFwd)$coefficients, digits = 4)
        cat("\n", "Alpha-to-enter = ", aEnter, ",   Alpha-to-remove = ", aRemove, "\n")
        if(plotRes == TRUE){
          plot(fitFwd, which=c(1,2))
        }
        return(invisible(fitFwd))
      } else {
        fitFwd <- fitBwd
      }
    }
  }
  else if(method=='stepwise'){
    # Step0 initiated
    fitBwd = stepfwd(fitBwd, fullmodel, forcedOut = forcedOut, aEnter = aEnter)
    step = 1
    num = length(attr(fullmodel$terms, 'dataClasses'))                               # maximum number of step
    while(-num < 0){
      print(paste0('Forward Selection Step: ', step))
      fitFwd = stepfwd(fitBwd, fullmodel, forcedOut = forcedOut, aEnter = aEnter)                     #function that tries to add a predictor to current model
      summary(fitFwd)
      if(identical(fitFwd, fitBwd) == T) {                                           #if no new predictors were added, it will stop
        cat("======== Final model ========", "\n")
        print(fitFwd$call)
        cat("Predictors forced in: ", forcedIn, "\n")
        cat("Predictors forced out: ", forcedOut, "\n", "\n")
        print(summary(fitFwd)$coefficients, digits = 4)
        cat("\n", "Alpha-to-enter = ", aEnter, ",  Alpha-to-remove = ", aRemove, "\n")
        return(invisible(fitFwd))
        break # no more addition
      }else {                                                                        #function that tries to remove a predictor from current model
        print(paste0('Backward Selection Step: ', step))
        fitBwd = stepbwd(fitFwd, fullmodel, forcedIn = forcedIn, aRemove = aRemove) 
        if(identical(fitFwd,fitBwd) != T){
          break # when there was a removal
        }
      }
      num = num + step
      step = step + 1
      }
  }
}

