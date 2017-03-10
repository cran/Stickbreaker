#' Fit all models to data
#'
#' Fit all models using \code{fit.models}. Used as a wrapper to combine outputs from \code{\link{estimate.d.MLE}}, \code{\link{estimate.d.RDB}}, \code{\link{estimate.d.sequential}}, \code{\link{fit.stick.model.given.d}}, \code{\link{fit.mult.model}}, \code{\link{fit.add.model}}, and \code{\link{summarize.fits.for.posterior.calc}}.
#'
#' @param data (Required). A list containing fitness values for allele combinations.
#' See examples in \code{data} or in the vignette.
#' @param d.range (Required). Range of realistic fitness changes from WT. "d" is the
#' difference between maximal fitness and WT fitness. For example, WT fitness for the
#' virus phix174 is about 20, but values up to about 25 have been observed, so d.range was set to:
#' d.range <- c(0.1,10).
#' @param d.adj.max (Optional). Default 1.1. When the estimation of d fails in
#' \code{\link{estimate.d.MLE}}, \code{\link{estimate.d.RDB}}, or \code{\link{estimate.d.sequential}} then d becomes 1.1*WT fitness so that models can be fit.
#' @param wts (Optional). Default c(2,1).
#'
#' @return List:\cr
#'  [[1]] \code{fit.smry} contains the R^2 and P values for all models
#'   coefficients \cr
#'  [[2]] \code{fit.stick} contains stickbreaking model details like the intercept, slope, and predicted fitness values for the given allele combinations.\cr
#'  [[2]] \code{fit.mult} contains multiplicative model details like the intercept, slope, and predicted fitness values for the given allele combinations.\cr
#'  [[2]] \code{fit.stick} contains additive model details like the intercept, slope, and predicted fitness values for the given allele combinations.\cr
#' @details Function is a wrapper that fits the data to all three models. First, it estimates $d$ and, using this estimate, fits the stickbreaking model. Then it fits the multiplicative and additive models.
#' @examples
#' fit.models(Chou.data, d.range=c(0.1, 10), d.adj.max=1.1, wts=c(2,1))
#' @seealso
#'  \code{\link{estimate.d.MLE}}, \code{\link{estimate.d.RDB}}, \code{\link{estimate.d.sequential}}, \code{\link{fit.stick.model.given.d}}, \code{\link{fit.mult.model}}, \code{\link{fit.add.model}}, \code{\link{summarize.fits.for.posterior.calc}}, \code{\link{estimate.d.MLE}}, \code{\link{estimate.d.RDB}}, \code{\link{estimate.d.sequential}}
#' @export
#'


fit.models <- function(data, d.range, d.adj.max=1.1, wts=c(2,1)){
  n.genos <- length(data[,1])
  n.muts <- length(data[1,])-1
  geno.matrix <- data[,seq(1, n.muts)]
  fit.matrix <- as.matrix(data[,(n.muts+1)])
  d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
  d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
  d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range)
  fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, run.regression=TRUE)
  fit.mult <- fit.mult.model(geno.matrix, fit.matrix, wts=wts)
  fit.add <- fit.add.model(geno.matrix, fit.matrix, wts=wts)
  fit.smry <- summarize.fits.for.posterior.calc(fit.stick, fit.mult, fit.add)
  return(list(fit.smry=fit.smry, fit.stick=fit.stick, fit.mult=fit.mult, fit.add=fit.add))
}

#' Fit the stickbreaking model to data for a given value of d
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}} or read in
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.stick.data}} or read in
#' @param d.here The value of d estimates are based on
#' @param wts Vector of weights to weight genotypes by. Used when
#'   \code{\link{generate.geno.weight.matrix}} is called (see that function).
#'   Default is \code{c(2,1)}, meaning weight single-mutation genotypes twice as heavily as others.
#'   Alternatively, vector of weights corresponding to geno.matrix can be provided.
#' @param run.regression \code{TRUE/FALSE} Run regression analysis when fitting model. See details.
#' @return List:\cr
#'  [[1]] \code{u.hats} are the estimated stickbreaking
#'   coefficients \cr
#'  [[2]] \code{R2} is proportion of fitness variation
#'   explained by model. Does not include wild type in calculation.\cr
#'  [[3]] \code{sig.hat} is estimate of sigma \cr
#'  [[4]] \code{logLike} is log-likelihood of the data under the fitted model. \cr
#'  [[5]] \code{regression.results} List of results when regressing effects of mutations against the background fitness
#' of mutations (see details). [[1]] \code{p.vals} gives p-value of each mutation, [[2]] \code{lm.intercepts} gives
#' estimated intercept for mutation, [[3]] \code{lm.slopes} gives slope for each mutation, [[4]] \code{P} is the
#' sum of the log of p-values. This is the summary statistic. [[5]] \code{fitness.of.backs} Matrix with fitness of backgrounds when each mutation (columns) is added to each genotype (rows).
#' [[6]] \code{effects.matrix} Matrix with fitness effect when given mutation (column) is added to given create genotype (row).
#' @details Note that the coefficient estimates are obtained by weighting. The
#'   default is to give wild type to single mutation genotypes twice the weight
#'   as all other comparisons based on the assumption that wild type is know
#'   with much lower error than the other genotypes. Alternatively, a vector of
#'   weights can be used with length the same as the number of genotypes in geno.matrix. \cr
#'
#'   In addition to R-squared we assess
#'   model fit by doing linear regression of background fitness against effect. When the model
#'   generating data and analyzing data are the same, the expected slope is zero and the p-values
#'   are uniform(0,1). The results from those regressions are returned in \code{regression.results}. \cr
#'  \code{run.regression} If you are doing simulations to assess parameter estimation only, you don't need to run
#'  regression. If you are using this function to generate data for model fitting, then this should be set to \code{TRUE}.
#'  @examples
#'  n.muts <- length(Khan.data[1,])-1
#'  geno.matrix <- Khan.data[,seq(1, n.muts)]
#'  fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
#'  d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))
#'  d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix,-100)$d.hat.RDB
#'  d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, c(0.1, 10), 1.1)
#'  fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, run.regression=TRUE)
#' @export


fit.stick.model.given.d <- function(geno.matrix, fit.matrix, d.here, wts=c(2,1), run.regression){

  # ---- Estimate coefficients ----
  n.genos <- dim(geno.matrix)[1]
  n.muts <- dim(geno.matrix)[2]
  weights.matrix <- generate.geno.weight.matrix(geno.matrix, fit.matrix, wts)
  fit.matrix.X.b <- rep(NA, dim(geno.matrix)[1])
  nonlog.vals <- 1-((fit.matrix[,1]-fit.matrix[1,1])/d.here)
  g.zero <- which(nonlog.vals >0)
  fit.matrix.X.b[g.zero] <- log(nonlog.vals[g.zero])
  #fit.matrix.X.b <- log(1-((fit.matrix[,1]-fit.matrix[1,1])/d.here))
  mut.all.geno.stick.X.b.coes <- matrix(nrow=length(geno.matrix[,1]), ncol=length(geno.matrix[1,]))
  mut.all.geno.stick.coes <- matrix(nrow=length(geno.matrix[,1]), ncol=length(geno.matrix[1,]))
  n.muts <- dim(geno.matrix)[2]
  geno.sim.strings.h <- apply(geno.matrix, MARGIN=1, FUN=paste, collapse="")
  u.hats.here <- rep(NA, n.muts)
  u.hats.here2 <- rep(NA, n.muts)
  names(u.hats.here) <- colnames(geno.matrix)
  #X.bs <- rep(NA, n.muts)
  mut.i.genos <- apply(geno.matrix, 2, function(x) which(x==1))
  if (is.list(mut.i.genos)==FALSE){
    mut.i.genos <- as.list(as.data.frame(mut.i.genos))
  }
  matrix.of.backgrounds <- matrix(data=NA, nrow=n.genos, ncol=n.muts)

  for (mut.i in 1:n.muts){
    for (geno.i in 1:length(mut.i.genos[[mut.i]])){   # so geno.i is indexing mut.i.genos (not geno.matrix or fit.matrix)
      #geno.ii <- mut.i.genos[geno.i, mut.i]
      geno.ii <- mut.i.genos[[mut.i]][geno.i]
      geno <- geno.matrix[geno.ii,]
      geno.background <- geno
      geno.background[mut.i] <- 0
      geno.back.string <- paste(geno.background, collapse="")
      back.id <- which(geno.sim.strings.h==geno.back.string)
      if (length(back.id)==1){
        matrix.of.backgrounds[geno.ii, mut.i] <- back.id
        mut.all.geno.stick.X.b.coes[mut.i.genos[[mut.i]][geno.i], mut.i] <- fit.matrix.X.b[geno.ii] - fit.matrix.X.b[back.id]
        mut.all.geno.stick.coes[mut.i.genos[[mut.i]][geno.i], mut.i] <- (fit.matrix[geno.ii,1] - fit.matrix[back.id,1])/(d.here - (fit.matrix[back.id,1]-fit.matrix[1,1]))
      }
    }
    non.na <- which(is.na(mut.all.geno.stick.X.b.coes[,mut.i])==FALSE)
    reweight <- weights.matrix[non.na,mut.i]/sum(weights.matrix[non.na, mut.i])
    u.hats.here[mut.i] <- 1-exp(sum(mut.all.geno.stick.X.b.coes[non.na,mut.i] * reweight))
    u.hats.here2[mut.i] <- sum(mut.all.geno.stick.coes[non.na,mut.i] * reweight)
    #X.bs[mut.i] <- sum(mut.all.geno.stick.X.b.coes[,mut.i] * weights.matrix[,mut.i], na.rm=TRUE)
  }

  # --- Remove infinite values from mut.all.geno.stick.coes---
  mut.all.geno.stick.coes[which(is.infinite(mut.all.geno.stick.coes), arr.ind=TRUE)] <- NA

  # --- get fitness of backgrounds ---
  fitness.of.backs <- matrix(nrow=n.genos, ncol=n.muts)
  for (mut.i in 1:n.muts){
    not.na <- which(is.na(matrix.of.backgrounds[,mut.i])==FALSE)
    not.na.backs <- matrix.of.backgrounds[not.na, mut.i]
    fitness.of.backs[not.na, mut.i] <- fit.matrix[not.na.backs]
  }

  # ---- Estimate sigma and assess model fit ----
  W.errors <- rep(NA, n.genos)
  W.g.pred <- apply(geno.matrix, 1, function(x) fit.matrix[1] +d.here*(1-prod(1-u.hats.here[which(x==1)])))
  W.errors <- fit.matrix - W.g.pred
  pred.matrix <- cbind(geno.matrix, string=geno.sim.strings.h, fit=fit.matrix, pred=W.g.pred, error=W.errors)

  # --- R2, sigma, lnL ---
  mean.fit <- mean(fit.matrix[2:n.genos,1], na.rm=TRUE)
  SS.null <- sum((fit.matrix[2:n.genos,1]-mean.fit)^2, na.rm=TRUE)
  SS.model <- sum(W.errors[2:n.genos]^2, na.rm=TRUE)
  R2.model <- 1-SS.model/SS.null
  sig.hat.nonlog <- sqrt(sum(W.errors[2:n.genos]^2, na.rm=T)/(length(which(is.na(fit.matrix)==FALSE))-2))
  lnL.nonlog <- sum(dnorm(W.errors, mean=0, sig.hat.nonlog, log=TRUE), na.rm=TRUE)

  # --- Regress background vs effect
  if (run.regression == TRUE){
    regression.results <- regress.back.fitness.vs.effect(fitness.of.backs, effects.matrix=mut.all.geno.stick.coes, n.muts=n.muts)
  } else{
    regression.results <- list(P=NA)
  }


  #results <- list(u.hats=u.hats.here, X.bs=X.bs, mut.all.geno.stick.X.b.coes=mut.all.geno.stick.X.b.coes, R2=R2.model, sig.hat=sig.hat.nonlog, Like=L.nonlog)
  results <- list(u.hats=u.hats.here, R2=R2.model, sig.hat=sig.hat.nonlog, logLike=lnL.nonlog, regression.results=regression.results, pred.matrix=pred.matrix)
  return(results)
}




#' Fit the multiplicative model to data
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}} or read in
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.mult.data}} or read in
#' @param wts Vector of weights to weight genotypes by. Used when
#'   \code{\link{generate.geno.weight.matrix}} is called (see that function).
#'   Default is \code{c(2,1)}.
#' @return List:\cr
#'  [[1]] \code{s.hats} are the estimated selection
#'   coefficients;\cr
#'  [[2]] \code{R2} is proportion of fitness variation
#'   explained by model. Does not include wild type in calculation.\cr
#'  [[3]] \code{sig.hat} is estimate of sigma \cr
#'  [[4]] \code{logLike} is log-likelihood of the data under the fitted model. \cr
#'  [[5]] \code{regression.results} List of results when regressing effects of mutations against the background fitness
#' of mutations (see details). [[1]] \code{p.vals} gives p-value of each mutation, [[2]] \code{lm.intercepts} gives
#' estimated intercept for mutation, [[3]] \code{lm.slopes} gives slope for each mutation, [[4]] \code{P} is the
#' sum of the log of p-values. This is the summary statistic. [[5]] \code{fitness.of.backs} Matrix with fitness of backgrounds when each mutation (columns) is added to each genotype (rows).
#' [[6]] \code{effects.matrix} Matrix with fitness effect when given mutation (column) is added to given create genotype (row).
#' @details \code{wts}:  The coefficient estimates are obtained by weighted comparisons. The
#'   default is to give wild type to single mutation genotype comparisons twice the weight
#'   as all other comparisons based on the assumption that wild type is know
#'   with much lower error than the other genotypes (actually it is assumed to be known with no error).
#'   @seealso \code{\link{fit.stick.model.given.d}}
#'   @examples
#'   fit.mult.model(
#'    data[,seq(1, length(Khan.data[1,])-1)],
#'    as.matrix(data[,(length(Khan.data[1,]))]),
#'    c(2,1))
#' @export


fit.mult.model <- function(geno.matrix, fit.matrix, wts=c(2,1)){

  w.wt <- fit.matrix[1,1]

  # ---- Estimate coefficients ----
  n.genos <- dim(geno.matrix)[1]
  n.muts <- dim(geno.matrix)[2]
  weights.matrix <- generate.geno.weight.matrix(geno.matrix, fit.matrix, wts)
  mut.all.geno.sel.coes.log <- matrix(nrow=length(geno.matrix[,1]), ncol=length(geno.matrix[1,]))
  mut.all.geno.sel.coes <- matrix(nrow=length(geno.matrix[,1]), ncol=length(geno.matrix[1,]))
  n.muts <- dim(geno.matrix)[2]
  geno.sim.strings.h <- apply(geno.matrix, MARGIN=1, FUN=paste, collapse="")
  s.hats.here <- rep(NA, n.muts)
  names(s.hats.here) <- colnames(geno.matrix)
  mut.i.genos <- apply(geno.matrix, 2, function(x) which(x==1))
  if (is.list(mut.i.genos)==FALSE){
    mut.i.genos <- as.list(as.data.frame(mut.i.genos))
  }
  matrix.of.backgrounds <- matrix(data=NA, nrow=n.genos, ncol=n.muts)

  for (mut.i in 1:n.muts){
    for (geno.i in 1:length(mut.i.genos[[mut.i]])){   # so geno.i is indexing mut.i.genos (not geno.matrix or fit.matrix)
      #geno.ii <- mut.i.genos[geno.i, mut.i]
      geno.ii <- mut.i.genos[[mut.i]][geno.i]
      geno <- geno.matrix[geno.ii,]
      geno.background <- geno
      geno.background[mut.i] <- 0
      geno.back.string <- paste(geno.background, collapse="")
      back.id <- which(geno.sim.strings.h==geno.back.string)
      if (length(back.id)==1){
        matrix.of.backgrounds[geno.ii, mut.i] <- back.id
        mut.all.geno.sel.coes[mut.i.genos[[mut.i]][geno.i], mut.i] <- (fit.matrix[geno.ii,1] - fit.matrix[back.id,1])/fit.matrix[back.id,1]
        mut.all.geno.sel.coes.log[mut.i.genos[[mut.i]][geno.i], mut.i]  <- log(fit.matrix[geno.ii, 1]) - log(fit.matrix[back.id, 1])  # this are observations on log scale
      }
    }
    non.na <- which(is.na(mut.all.geno.sel.coes.log[,mut.i])==FALSE)
    reweight <- weights.matrix[non.na,mut.i]/sum(weights.matrix[non.na, mut.i])
    s.hats.here[mut.i] <- exp(sum(mut.all.geno.sel.coes.log[non.na, mut.i]*reweight)) - 1
  }

  # --- get fitness of backgrounds ---
  fitness.of.backs <- matrix(nrow=n.genos, ncol=n.muts)
  for (mut.i in 1:n.muts){
    not.na <- which(is.na(matrix.of.backgrounds[,mut.i])==FALSE)
    not.na.backs <- matrix.of.backgrounds[not.na, mut.i]
    fitness.of.backs[not.na, mut.i] <- fit.matrix[not.na.backs]
  }

  # ---- Estimate sigma and assess model fit ----
  W.errors <- rep(NA, n.genos)
  geno.coes <- t(apply(geno.matrix, 1, function(x) x*s.hats.here+1))
  W.g.pred <- as.matrix(apply(geno.coes, 1, function(x) w.wt*prod(x)))
  W.errors <- fit.matrix - W.g.pred
  pred.matrix <- cbind(geno.matrix, string=geno.sim.strings.h, fit=fit.matrix, pred=W.g.pred, error=W.errors)

  # --- R2, sigma, lnL ---
  mean.fit <- mean(fit.matrix[2:n.genos,1], na.rm=TRUE)
  SS.null <- sum((fit.matrix[2:n.genos,1]-mean.fit)^2, na.rm=TRUE)
  SS.model <- sum(W.errors[2:n.genos]^2, na.rm=TRUE)
  R2.model <- 1-SS.model/SS.null
  sig.hat.nonlog <- sqrt(sum(W.errors[2:n.genos]^2, na.rm=T)/(length(which(is.na(fit.matrix)==FALSE))-2))
  lnL.nonlog <- sum(dnorm(W.errors, mean=0, sig.hat.nonlog, log=TRUE), na.rm=TRUE)

  # --- Regress background vs effect
  regression.results <- regress.back.fitness.vs.effect(fitness.of.backs, effects.matrix=mut.all.geno.sel.coes, n.muts=n.muts)

  results <- list(s.hats=s.hats.here, R2=R2.model, sig.hat=sig.hat.nonlog, logLike=lnL.nonlog, regression.results=regression.results, pred.matrix=pred.matrix)
  return(results)
}






#' Fit the additive model to data
#'
#' @param geno.matrix Genotype matrix generated in
#'   \code{\link{generate.geno.matrix}} or read in
#' @param fit.matrix Fitness matrix generated in
#'   \code{\link{sim.mult.data}} or read in
#' @param wts Vector of weights to weight genotypes by. Used when
#'   \code{\link{generate.geno.weight.matrix}} is called (see that function).
#'   Default is \code{c(2,1)}.
#' @return List:\cr
#'   [[1]] \code{w.hats} are the estimated additive effect
#'   coefficients;\cr
#'   [[2]] \code{R2} is proportion of fitness variation
#'   explained by model. Does not include wild type in calculation.\cr
#'   [[3]] \code{sig.hat} is estimate of sigma \cr
#'   [[4]] \code{logLike} is log-likelihood of the data under the fitted model. \cr
#'  [[5]] \code{regression.results} List of results when regressing effects of mutations against the background fitness
#' of mutations (see details). [[1]] \code{p.vals} gives p-value of each mutation, [[2]] \code{lm.intercepts} gives
#' estimated intercept for mutation, [[3]] \code{lm.slopes} gives slope for each mutation, [[4]] \code{P} is the
#' sum of the log of p-values. This is the summary statistic. [[5]] \code{fitness.of.backs} Matrix with fitness of backgrounds when each mutation (columns) is added to each genotype (rows).
#' [[6]] \code{effects.matrix} Matrix with fitness effect when given mutation (column) is added to given create genotype (row).
#' @details \code{wts}:  The coefficient estimates are obtained by weighted comparisons. The
#'   default is to give wild type to single mutation genotype comparisons twice the weight
#'   as all other comparisons based on the assumption that wild type is know
#'   with much lower error than the other genotypes (actually it is assumed to be known with no error).
#'   @examples
#'   n.muts <- length(Khan.data[1,])-1
#'   geno.matrix <- Khan.data[,seq(1, n.muts)]
#'   fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
#'   fit.add.model(geno.matrix, fit.matrix, c(2,1))
#'   @seealso \code{\link{fit.stick.model.given.d}}
#' @export

fit.add.model <- function(geno.matrix, fit.matrix, wts=c(2,1)){

  w.wt <- fit.matrix[1,1]
  # ---- Estimate coefficients ----
  n.genos <- dim(geno.matrix)[1]
  n.muts <- dim(geno.matrix)[2]
  weights.matrix <- generate.geno.weight.matrix(geno.matrix, fit.matrix, wts)
  mut.all.geno.add.coes <- matrix(nrow=length(geno.matrix[,1]), ncol=length(geno.matrix[1,]))
  n.muts <- dim(geno.matrix)[2]
  geno.sim.strings.h <- apply(geno.matrix, MARGIN=1, FUN=paste, collapse="")
  w.hats.here <- rep(NA, n.muts)
  names(w.hats.here) <- colnames(geno.matrix)
  mut.i.genos <- apply(geno.matrix, 2, function(x) which(x==1))
  if (is.list(mut.i.genos)==FALSE){
    mut.i.genos <- as.list(as.data.frame(mut.i.genos))
  }
  matrix.of.backgrounds <- matrix(data=NA, nrow=n.genos, ncol=n.muts)

  for (mut.i in 1:n.muts){
    for (geno.i in 1:length(mut.i.genos[[mut.i]])){   # so geno.i is indexing mut.i.genos (not geno.matrix or fit.matrix)
      #geno.ii <- mut.i.genos[geno.i, mut.i]
      geno.ii <- mut.i.genos[[mut.i]][geno.i]
      geno <- geno.matrix[geno.ii,]
      geno.background <- geno
      geno.background[mut.i] <- 0
      geno.back.string <- paste(geno.background, collapse="")
      back.id <- which(geno.sim.strings.h==geno.back.string)
      if (length(back.id)==1){
        matrix.of.backgrounds[geno.ii, mut.i] <- back.id
        mut.all.geno.add.coes[mut.i.genos[[mut.i]][geno.i], mut.i]  <- (fit.matrix[geno.ii, 1]) - (fit.matrix[back.id, 1])  # this are observed effects
      }
    }
    non.na <- which(is.na(mut.all.geno.add.coes[,mut.i])==FALSE)
    reweight <- weights.matrix[non.na,mut.i]/sum(weights.matrix[non.na, mut.i])
    w.hats.here[mut.i] <- sum(mut.all.geno.add.coes[non.na, mut.i]*reweight)
  }

  # --- get fitness of backgrounds ---
  fitness.of.backs <- matrix(nrow=n.genos, ncol=n.muts)
  for (mut.i in 1:n.muts){
    not.na <- which(is.na(matrix.of.backgrounds[,mut.i])==FALSE)
    not.na.backs <- matrix.of.backgrounds[not.na, mut.i]
    fitness.of.backs[not.na, mut.i] <- fit.matrix[not.na.backs]
  }

  # ---- Estimate sigma and assess model fit ----
  W.errors <- rep(NA, n.genos)
  geno.coes <- t(apply(geno.matrix, 1, function(x) x*w.hats.here))
  w.wt <- fit.matrix[1,1]
  W.g.pred <- as.matrix(apply(geno.coes, 1, function(x) w.wt+sum(x)))
  W.errors <- fit.matrix - W.g.pred

  pred.matrix <- cbind(geno.matrix, string=geno.sim.strings.h, fit=fit.matrix, pred=W.g.pred, error=W.errors)

  # --- R2, sigma, lnL ---
  mean.fit <- mean(fit.matrix[2:n.genos,1], na.rm=TRUE)
  SS.null <- sum((fit.matrix[2:n.genos,1]-mean.fit)^2, na.rm=TRUE)
  SS.model <- sum(W.errors[2:n.genos]^2, na.rm=TRUE)
  R2.model <- 1-SS.model/SS.null
  sig.hat.nonlog <- sqrt(sum(W.errors[2:n.genos]^2, na.rm=T)/(length(which(is.na(fit.matrix)==FALSE))-2))
  lnL.nonlog <- sum(dnorm(W.errors, mean=0, sig.hat.nonlog, log=TRUE), na.rm=TRUE)

  # --- Regress background vs effect
  regression.results <- regress.back.fitness.vs.effect(fitness.of.backs, effects.matrix=mut.all.geno.add.coes, n.muts=n.muts)

  results <- list(w.hats=w.hats.here, R2=R2.model, sig.hat=sig.hat.nonlog, logLike=lnL.nonlog, regression.results=regression.results, pred.matrix=pred.matrix)
  return(results)
}



#' Linear regression of background fitness against effects
#'
#' @param fitness.of.backs Matrix with background fitness for \code{effects.matrix}
#' @param effects.matrix Effect of mutation in that genotype
#' @param n.muts Number of mutations
#' @return List: \cr
#' [[1]] \code{p.vals} P-value for linear regression for each mutation \cr
#' [[2]] \code{lm.intercepts} Intercepts for each mutation
#' [[3]] \code{lm.slopes} Slopes for each mutation
#' [[4]] \code{P} Sum of the logs of the p-values.
#' [[5]] \code{fitness.of.backs} Matrix with fitness of backgrounds when each mutation (columns) is added to each genotype (rows).
#' [[6]] \code{effects.matrix} Matrix with fitness effect when given mutation (column) is added to given create genotype (row).
#' If there is insufficient data to do regression, a warning is returned.
#' @details For each mutation, function does simple linear regression using \code{lm()}.
#' The product of logged p-values (\code{P}) is the summary statistic used in model selection.
#' When the model that generated the data and the model analyzing the data match, the expected slope
#' of the regression line is zero.
#' @export


regress.back.fitness.vs.effect <- function(fitness.of.backs, effects.matrix, n.muts){
  p.vals <- rep(NA, n.muts)
  lm.intercepts <- rep(NA, n.muts)
  lm.slopes <- rep(NA, n.muts)

  for (mut.i in 1:n.muts){
    not.na <- which(is.na(fitness.of.backs[,mut.i])==FALSE)
    b.fits <- fitness.of.backs[not.na, mut.i]
    effects <- effects.matrix[not.na, mut.i]
    s.lm <- summary(lm(effects ~ b.fits))

    if (is.null(s.lm$fstatistic[1])==FALSE){
      p.vals[mut.i] <- pf(s.lm$fstatistic[1], s.lm$fstatistic[2], s.lm$fstatistic[3], lower.tail=FALSE)
      lm.intercepts[mut.i] <- s.lm$coefficients[1,1]
      lm.slopes[mut.i] <- s.lm$coefficients[2,1]
      #p.smry <- sum(log(p.vals),na.rm=TRUE)
      #p.smry.2 <- prod(1-p.vals)
    } else{
      p.vals[mut.i] <- NA
      lm.intercepts[mut.i] <- NA
      lm.slopes[mut.i] <- NA
    }
  }
  p.smry <- sum(log(p.vals),na.rm=TRUE)
  p.smry.2 <- prod(1-p.vals)
  return(list(p.vals=p.vals, lm.intercepts = lm.intercepts, lm.slopes=lm.slopes, P=p.smry, fitness.of.backs=fitness.of.backs, effects.matrix=effects.matrix))
}


#' Extracts summary statistics from each model needed for posterior calculation
#'
#' @param fit.stick List returned from fitting data to stickbreaking model in \code{fit.stick.model.given.d}
#' @param fit.mult List returned from fitting data to multiplicative model in \code{fit.mult.model}
#' @param fit.add List returned from fitting data to additive model in \code{fit.add.model}
#' @return Vector of summary statistics under each model
#' @details Extracts R2, sigma, log-likelihood and P-statistic under stickbreaking, multiplicative and additive models
#' @export

summarize.fits.for.posterior.calc <- function(fit.stick, fit.mult, fit.add){
  f.stick <- c(R2.stick=fit.stick$R2, P.stick=round(fit.stick$regression.results$P,5))
  f.mult <- c(R2.mult=fit.mult$R2, P.mult=round(fit.mult$regression.results$P,5))
  f.add <- c(R2.add=fit.add$R2, P.add=round(fit.add$regression.results$P,5))

  fit.smry <- c(f.stick, f.mult, f.add)
  fit.smry <- as.data.frame(t(as.matrix(fit.smry)))
  return(fit.smry)
}
