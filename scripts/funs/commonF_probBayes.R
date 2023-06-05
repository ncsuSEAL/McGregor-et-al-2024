calcProb <- function(z, prior=NULL){
  like <- densEDist(z)
  like0 <- densEStable(z)
  
  # if(is.na(like)) like <- 0
  # if(is.na(like0)) like0 <- 0
  
  if(is.null(prior)){
    a <- like
    b <- like0
  } else {
    a <- like * prior
    b <- like0 * (1-prior) 
  }
  return(a / (a+b))
}

applyBayes <- function(allProbs, pt=NULL, binPt, runType="train"){
  if(runType=="app"){
    prior <- allProbs
  } else if(runType=="train"){
    prior <- allProbs[pointid==pt, prior]
  }
  
  probsBayes <- sapply(binPt$ewma, calcProb, prior)
  
  ## have to do a correction here: some bad model fits (e.g. S2 in chatthin)
  ## lead to very high + ewma values. These are outside the bounds of the
  ## density function, so NA is returned. Essentially this is 0 likelihood,
  ## but for mathematical consistency we're setting it as just under the final
  ## bound for the function (e.g. densF(5.2) = 3.2e-6)
  probsBayes[is.na(probsBayes)] <- 1e-6
  return(probsBayes)
}



# priorSim <- function(){
#   set.seed(456)
#   # zewmaSim <- binPt$ewma
#   zewmaSim <- c(runif(80), runif(20, min=-3, max=-2))
#   
#   priorSeq <- c(0.01, 0.5, 0.99)
#   lty <- c(1,2,3)
#   col <- c("black")
#   
#   ## non-Bayes
#   plot(NA, xlim=c(1, length(zewmaSim)), ylim=c(0,1), main="Non-Bayes",
#        ylab="Probability", xlab="")
#   bleh <- sapply(1:length(priorSeq), function(p){
#     probsBayes <- sapply(zewmaSim, function(z){
#       return(calcProb(z, prior=NULL))
#     })
#     lines(probsBayes, col=col[1], lty=lty[p])
#     return(probsBayes)
#   })
#   
#   ## Bayes
#   plot(NA, xlim=c(1, length(zewmaSim)), ylim=c(0,1), main="Bayesian",
#        ylab="Probability", xlab="")
#   legend("topleft", legend=paste0("Prior=", priorSeq), lty=lty, col="black",
#          bty="n")
#   bleh <- sapply(1:length(priorSeq), function(p){
#     probsBayes <- sapply(zewmaSim, function(z){
#       return(calcProb(z, priorSeq[p]))
#     })
#     lines(probsBayes, col=col[1], lty=lty[p])
#     return(probsBayes)
#   })
# }
# 
# layout(matrix(1:2, ncol=2))
# priorSim()

