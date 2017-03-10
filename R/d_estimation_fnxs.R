#' Wrapper function so log-likelihood of stickbreaking can be extracted by optimize() function
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}}
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.stick.data}}
#' @param d.vect Parameter to be maximized over
#' @param wts Vector of weights to weight genotypes by. Used when
#'   \code{\link{generate.geno.weight.matrix}} is called (see that function).
#'   Default is \code{c(2,1)}, meaning weight single-mutation genotypes twice as heavily as others.
#'   Alternatively, vector of weights corresponding to geno.matrix can be provided.
#' @return Vector of log-likelihoods
#' @details Calls \code{fit.stick.model.given.d} to get log likelihood
#' @export

calc.stick.logLn <- function(geno.matrix, fit.matrix, d.vect, wts=c(2,1)){
  L.vect <- sapply(d.vect, function(x) fit.stick.model.given.d(geno.matrix, fit.matrix, x, wts, run.regression=FALSE)$logLike)
  return(L.vect)
}


#' Find MLE of d
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}}
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.stick.data}}
#' @param d.range Interval of d to search for maximum over
#' @param accuracy \code{tol} to send \code{optimize} function
#' @param wts Vector of weights to weight genotypes by. Used when
#'   \code{\link{generate.geno.weight.matrix}} is called (see that function).
#'   Default is \code{c(2,1)}, meaning weight single-mutation genotypes twice as heavily as others.
#'   Alternatively, vector of weights corresponding to geno.matrix can be provided.
#' @return MLE of d
#' @details Maximizes the function \code{\link{calc.stick.logLn}} using \code{\link{optimize}}
#' @examples
#' n.muts <- length(Khan.data[1,])-1
#' geno.matrix <- Khan.data[,seq(1, n.muts)]
#' fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
#' estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))

#' @export

estimate.d.MLE <- function(geno.matrix, fit.matrix, d.range, accuracy=0.001, wts=c(2,1)){
  d.vect <- seq(d.range[1], d.range[2], length.out=100)
  maxLike <- optimize(f=calc.stick.logLn, interval=d.range, geno.matrix=geno.matrix, fit.matrix=fit.matrix, wts=wts, maximum=TRUE, tol=accuracy)
  return(maxLike$maximum)
}




#' Estimate d using relative distance to boundary (RDB) methods
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}} or read in.
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.stick.data}} or read in.
#' @param no.est What to return when no estimate is obtained. Default is -100. Throws warning unless a number of NA is used.
#' @return List: \cr
#'    [[1]] \code{d.hat} Is the best RDB estimate. It is the
#'   median of the positive individual RDB values. \cr
#'   [[2]] \code{d.hat.RDB.all} RDB estimate based on all individual estimates (not only postie values). \cr
#'   [[3]] \code{d.hat.RDB.ind} Vector of individual RDB values (see desalt).\cr
#'   [[4]] \code{d.hat.RDB.other} has alternative ways to combine the individual
#'   RDB estimates (i.e. median, mean, all values, postie values only).\cr
#'   \code{d.hat.median.pos} is preferred estimator (\code{d.hat} above) from
#'   the median of positive individual values. \cr \code{d.hat.median} is median
#'   of all values. \cr \code{d.hat.mean} is mean of all values.\cr
#'   \code{d.hat.mean.pos} is mean of positive values.
#' @details The method calculates RDB for each genotype and its complement. The
#'   \code{d.hat.RDB.all} indicate the genotype pair that produces the estimate.
#' @examples
#' n.muts <- length(Khan.data[1,])-1
#' geno.matrix <- Khan.data[,seq(1, n.muts)]
#' fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
#' estimate.d.RDB(geno.matrix, fit.matrix,-100)
#' @export


estimate.d.RDB <- function(geno.matrix, fit.matrix, no.est=-100){   # calculates the relative distance to boundary estimators
  n.genos <- dim(geno.matrix)[1]
  n.muts <- dim(geno.matrix)[2]
  X.b <- rep(NA, n.genos)
  W.b <- fit.matrix[2:(n.genos-1),1]
  W.wt <- fit.matrix[1,1]
  W.k <- fit.matrix[n.genos, 1]
  geno.strings <- rep(NA, n.genos)
  singles <- which(apply(geno.matrix, 1, function(x) sum(x))==1)
  geno.strings <- apply(geno.matrix, MARGIN=1, FUN=paste, collapse="")
  geno.comp <- abs(geno.matrix-1)  # get complement of genotype
  geno.comp.strings <- apply(geno.comp, 1, FUN=paste, collapse="")

  W.b.c <- rep(NA, n.genos)
  for (geno.i in 2:(n.genos-1)){
    W.b.c[geno.i] <- fit.matrix[match(geno.comp.strings[geno.i], geno.strings),1]
  }
  W.b.c <- W.b.c[-n.genos]
  W.b.c <- W.b.c[-1]
  e.X.b <- (W.k-W.b)/(W.b.c-W.wt)
  d.hat <- (W.b-W.wt)/(1-e.X.b)
  pos <- which(e.X.b > 0 & e.X.b < 1)

  if (length(pos)>0){
    d.hat.pos <- (W.b[pos]-W.wt)/(1-e.X.b[pos])
  } else{
    d.hat.pos <- no.est
  }

  geno.name.both <- cbind(geno.strings[2:(n.genos-1)], geno.comp.strings[2:(n.genos-1)])
  names(d.hat) <- apply(geno.name.both, 1, function(x) paste(x, collapse="."))

  d.hat.mean <- mean(d.hat, na.rm=TRUE)
  d.hat.median <- median(d.hat, na.rm=TRUE)
  d.hat.median.pos <- median(d.hat.pos, na.rm=TRUE)
  d.hat.mean.pos <- mean(d.hat.pos, na.rm=TRUE)
  d.hat.med.singles <- median(d.hat[singles], na.rm=T)
  d.hat.RDB.others <- list(d.hat.median.pos=d.hat.median.pos, d.hat.median=d.hat.median, d.hat.mean=d.hat.mean, d.hat.mean.pos=d.hat.mean.pos)
  results <- list(d.hat.RDB = d.hat.median.pos, d.hat.RDB.all=d.hat.median, d.hat.RDB.ind=d.hat, d.hat.RDB.others=d.hat.RDB.others)
  #results <- list(d.hat.median, d.hat.mean, d.hat.med.singles, d.hat, d.hat.mean.pos, d.hat.median.pos)
  #names(results) <- c("d.hat.median", "d.hat.mean", "d.hat.med.sing", "d.hat.by.geno", "d.hat.mean.pos", "d.hat.median.pos")
  return(results)
}




#' Estimate d using sequential method
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}} or read in.
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.stick.data}} or read in.
#' @param d.hat.MLE The MLE estimate of d from \code{\link{estimate.d.MLE}}
#' @param d.hat.RDB The estimate of d from \code{\link{estimate.d.RDB}}
#' @param d.range Interval of d where MLE was searched. Defines when valid MLE exists (see details).
#' @param d.max.adj When forced to use the maximum estimator, the estimate is adjusted upwards
#' by this factor (see details). Default = 1.1 (inflate observation 10\%).
#' @return Estimate of d. Name indicates methods used.
#' @details If a valid MLE exists, function returns it. When MLE is at the boundary
#' as defined by \code{d.range} or less than the observed distance from wild type fitness to the maximum fitness,
#' this is not considered a valid estimate. If MLE is not valid, but RDB
#' estimate exists, it returns RDB estimate. If neither exists, returns the Max estimated.
#' The Max estimate is based on the largest observed fitness times a factor \code{d.max.adj}.
#' Name of the returned object indicates the method estimate is based on (MLE, RDB or Max).
#' @examples
#' n.muts <- length(Khan.data[1,])-1
#' geno.matrix <- Khan.data[,seq(1, n.muts)]
#' fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
#' d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))
#' d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix,-100)$d.hat.RDB
#' estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, c(0.1, 10), 1.1)

#' @export

estimate.d.sequential <- function(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range, d.max.adj=1.1){
  d.to.max <- max(fit.matrix[,1], na.rm=TRUE) -fit.matrix[1,1]
  if (is.na(d.hat.MLE)==FALSE & d.hat.MLE < 0.99*d.range[2] & d.hat.MLE > 1.01*d.range[1] & d.hat.MLE > d.to.max){
    d.hat.seq <- d.hat.MLE
    names(d.hat.seq) <- "MLE"
  } else if (is.na(d.hat.RDB)==FALSE & d.hat.RDB < 0.99*d.range[2] & d.hat.RDB > 1.01*d.range[1] &  d.hat.RDB > d.to.max){
    d.hat.seq <- d.hat.RDB
    names(d.hat.seq) <- "RDB"
  } else{
    d.hat.seq <- max(fit.matrix, na.rm=TRUE)*d.max.adj - fit.matrix[1]
    names(d.hat.seq) <- "max"
  }
  return(d.hat.seq)
}
