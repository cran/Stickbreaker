#' Generate genotype matrix for given number of mutations.
#'
#' @param n Number of mutations.
#' @return Matrix with all genotypes made of \code{n} mutations. Each column is a mutation (0=absent, 1=present), rows are genotypes.
#' @details Used when simulating data for a complete network. Generally called internally during simulations to create matrix for which fitness values will then be simulated.
#' Minimum value of \code{n} is 2.
#' @examples generate.geno.matrix(5)
#' @seealso \code{\link{sim.stick.data}}
#' @export


generate.geno.matrix <- function(n){
  geno.matrix.sim <- matrix(nrow=2^n, ncol=n)
  g <- matrix(data=NA, nrow=4, ncol=2)
  g[1,] <- c(0,0)
  g[2,] <- c(1,0)
  g[3,] <- c(0,1)
  g[4,] <- c(1,1)

  if (n >= 3){
    for (n.added in 3:n){
      n.cur <- length(g[,1])
      g <- rbind(g,g)
      c.to.add <- c(rep(0, n.cur), rep(1, n.cur))
      g <- cbind(g, c.to.add)
    }
  }
  colnames(g) <- paste("mut",seq(1, n), sep=".")
  return(g)
}

#' Simulate data under stickbreaking model.
#'
#' @param n.muts Number of mutations.
#' @param coes Vector of stickbreaking coefficients for each mutation
#' @param sigma Noise to add to expectations
#' @param d.true The distance to the boundary to use in simulations (d)
#' @param w.wt Fitness of the wildtype
#' @param sim.geno.matrix If TRUE, simulate geno.matrix. If FALSE, use geno.matrix provided.
#' @param geno.matrix Specifies which genotypes to simulated.
#' @return List:\cr
#' [[1]] \code{fit.matrix} is matrix with (simulated) observed data \cr
#' [[2]] \code{fit.matrix.exp} is same matrix without error
#' @details Function takes coefficients (\code{coes}) and generates the expected
#' fitness values for each genotype. This set of expected values is given in \code{fit.matrix.exp}.
#' Then it adds normal error to them to produce the "observed" data (\code{fit.matrix}).
#' @examples
#'  geno.matrix <- generate.geno.matrix(5)
#'  coes <- rep(0.1, 5)
#'  sim.stick.data(5, coes, 0.1, 1, 1, geno.matrix)
#' @export


sim.stick.data <- function(n.muts, coes, sigma, d.true, w.wt, geno.matrix,sim.geno.matrix){
  n.genos <- length(geno.matrix[,1])
  geno.coes <- geno.matrix*coes
  fit.matrix.exp <- as.matrix(apply(geno.coes, 1, function(x) w.wt + d.true*(1-prod(1-x))))
  error <- c(0, rnorm(n=n.genos-1, mean=0, sd=sigma))
  fit.matrix <- as.matrix(fit.matrix.exp + error)
  results <- list(fit.matrix.exp=fit.matrix.exp, fit.matrix=fit.matrix)
  return(results)
}


#' Simulate and fit batch data under stickbreaking model
#'
#' @param mut.vals Vector of number of mutations to simulate
#' @param coe.vals Vector of stickbreaking coefficients to simulate
#' @param sig.vals Vector of sigma values to simulate
#' @param d.true The distance to the boundary to use in simulations (d)
#' @param d.range Vector of range values. When estimating d under MLE, what range should be searched (see details). Default is \code{c(0.1, 10)}.
#' @param w.wt Fitness of the wildtype
#' @param n.reps.ea Number of replicates per parametric condition
#' @param print.status \code{TRUE/FALSE}. Should loop counters be printed.
#' @param fit.methods Vector of all methods of estimating d to then fit model and output results.
#' Accepts "MLE", "RDB", "max", "seq", "RDB.all" and "All". "All" does all methods. Default is "seq". Case sensitive.
#' @param outdir The path to write output files to (see details about file names).
#' @param wts The weight assigned to wildtype vs other genotypes when estimating parameters (see details). Default c(2,1).
#' @param d.max.adj When forced to use the maximum estimator, the estimate is adjusted upwards
#' by this factor (see details). Default = 1.1 (inflate observation 10\%).
#' @param run.regression \code{TRUE/FALSE} Run regression analysis when fitting model. See details.
#' @param RDB.method Indicates which RDB method to use when doing sequential estimation. Options are "pos" and "all".
#' "pos" option is better when mutations are strictly beneficial: "all" is appropriate when some or all mutations are deleterious.
#' @details Function contains a loop for combining each value of mut.vals, coe.vals and sig.vals,
#' generating data under the stickbreaking model and then fitting it. The fit.methods argument allows
#' user to evaluate performance of multiple methods at one time. Note that estimate of d under MLE is
#' restricted to \code{d.range}. Using a reasonable upper bound here is valuable so that the stickbreaking model
#' remains distant from the additive model (i.e. as d gets large and the stickbreaking coefficients get small,
#' the stickbreaking model converges to the additive model). \cr\cr
#' Results are written to files; the name of the output files are formed by concatenating the outpath
#' argument to the item in the fit.methods. Separate files are generated for each method (e.g. MLE, RDB, seq).
#' Separate files are also generated for each
#' number of mutations (because the dimensionality of the output file changes with the number of mutations).
#' The output files are tab-delimited text files with one row per replicate.
#' The first 5 columns provide the parameter values and the rest of the columns give parameter estimates and
#' measures of fit. \cr\cr
#' \code{wts}:  The coefficient estimates are obtained by weighted comparisons. The
#' default is to give wild type to single mutation genotype comparisons twice the weight
#' as all other comparisons based on the assumption that wild type is know
#' with much lower error than the other genotypes (actually it is assumed to be known with no error). \cr\cr
#'  \code{run.regression} If you are doing simulations to assess parameter estimation only, you don't need to run
#'  regression. If you are using this function to generate data for model fitting, then this should be set to \code{TRUE}.
#' @return Nothing. Instead results are written to output files and deposited in inst/extdata.
#' The files are named by appending the method
#' @examples
#' \dontrun{
#' sim.fit.stick.data.batch(c(3,4,5),
#'   c(0.1, 0.3, 0.5),
#'   c(0.02, 0.05, 0.08),
#'   1,
#'   c(0.1, 10),
#'   1,
#'   10,
#'   print.status=FALSE,
#'   fit.methods="seq",
#'   outdir="~/Desktop",
#'   c(2,1),
#'   1.0,
#'   run.regression="FALSE",
#'   RDB.method="pos")
#'   }
#' @export
#'


sim.fit.stick.data.batch <- function(mut.vals, coe.vals, sig.vals, d.true, d.range=c(0.1, 10), w.wt, n.reps.ea=100, print.status=FALSE, fit.methods=c("seq"), outdir, wts=c(2,1), d.max.adj=1.1, run.regression=TRUE, RDB.method){
  outfile <- "Stick_batch_out_test"
  outpath <- paste(outdir, outfile, sep="/")
  if (fit.methods == "All"){
    fit.methods2 <- c("MLE", "RDB", "max", "seq", "RDB.all")
  } else{
    fit.methods2 <- fit.methods
  }
  for (mut.i in 1:length(mut.vals)){
    n.muts <- mut.vals[mut.i]
    geno.matrix <- generate.geno.matrix(n.muts)
    n.genos <- 2^n.muts
    if (print.status==TRUE) print(paste("n.muts =", n.muts))
    first.results <- rep(TRUE, length(fit.methods2))

    for (coe.i in 1:length(coe.vals)){
      coes <- rep(coe.vals[coe.i], n.muts)
      if (print.status ==TRUE) print(paste("coe = ", coes[1]))

      for (sig.i in 1:length(sig.vals)){
        sigma <- sig.vals[sig.i]
        if (print.status == TRUE) print(paste("sigma = ", sigma))

        for (rep.i in 1:n.reps.ea){
          stick.sim.data <- sim.stick.data(n.muts, coes, sigma, d.true, w.wt, geno.matrix)
          fit.matrix <- stick.sim.data$fit.matrix
          d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
          d.hat.RDB.list <- estimate.d.RDB(geno.matrix, fit.matrix)
          d.hat.RDB <- d.hat.RDB.list$d.hat.RDB
          d.hat.RDB.all <- d.hat.RDB.list$d.hat.RDB.all
          if (RDB.method == "pos"){
            RDB.to.send <- d.hat.RDB
          } else if (RDB.method == "all"){
            RDB.to.send <- d.hat.RDB.all
          }
          d.hat.max <- max(fit.matrix, na.rm=TRUE)*d.max.adj - w.wt
          d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, RDB.to.send, d.range)
          d.hat.list <- list(d.hat.MLE=d.hat.MLE, d.hat.RDB=d.hat.RDB, d.hat.max=d.hat.max, d.hat.seq=d.hat.seq, d.hat.RDB.all=d.hat.RDB.all)
          parm.row <- unlist(list(model="stick", n.muts=n.muts, coes=coes[1], sigma=sigma, d.true=d.true, rep=rep.i))

          for (method.i in 1:length(fit.methods2)){
            method.dex <- which(names(d.hat.list)==paste("d.hat.", fit.methods2[method.i], sep=""))
            fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.list[[method.dex]], wts=wts, run.regression)
            method <- fit.methods2[method.i]
            if (method=="seq"){
              method.dhat <- names(d.hat.seq)
            } else{
              method.dhat <- method  #method[method.dex]
            }
            fit.row <- unlist(list(method=method.dhat, d.hat=round(d.hat.list[[method.dex]],5), round(unlist(fit.stick[1:4]),5), P=round(fit.stick$regression.results$P,5)))
            one.row <- as.data.frame(t(c(parm.row, fit.row)))
            output.file <- paste(outpath, "_", method, "_", n.muts, "muts.txt", sep="")
            if (first.results[method.i] == TRUE){
              write.table(x=one.row, file=output.file, append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
              if (print.status ==TRUE)  print(paste("Writing results to: ", output.file), sep="")
              #con <- file(description=output.file)
              #open(con)
              first.results[method.i] <- FALSE
            } else{
              write.table(x=one.row, file=output.file, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
              #writeLines(text=one.row, con=con, sep="\t")
              #cat(one.row, file=con, sep="\t")
            }
          }
        }  #next rep.i
      } #next sig.i
    } # next coe.i
  } # next mut.i
}


#' Simulate data under multiplicative model.
#'
#' @param n.muts Number of mutations.
#' @param selcoes Vector of selection coefficients for each mutation
#' @param sigma Noise to add to expectations
#' @param w.wt Fitness of the wildtype
#' @param geno.matrix Specifies which genotypes are simulated.
#' @return List:\cr
#' [[1]] \code{fit.matrix} is matrix with (simulated) observed data \cr
#' [[2]] \code{fit.matrix.exp} is same matrix without error
#' @details Function takes selection coefficients (\code{selcoes}) and generates the expected
#' fitness values for each genotype. This set of expected values is returned as \code{fit.matrix.exp}.
#' Then it adds normal error to them to produce the "observed" data (\code{fit.matrix}).
#' @examples
#'  geno.matrix <- generate.geno.matrix(3)
#'  selcoes <- rep(0.3, 3)
#'  sim.mult.data(3, selcoes, 0.1, 1, geno.matrix)
#' @export


sim.mult.data <- function(n.muts, selcoes, sigma, w.wt, geno.matrix){
  n.genos <- length(geno.matrix[,1])
  geno.coes <- geno.matrix*selcoes +1
  fit.matrix.exp <- as.matrix(apply(geno.coes, 1, function(x) w.wt*prod(x)))
  error <- c(0, rnorm(n=n.genos-1, mean=0, sd=sigma))
  fit.matrix <- as.matrix(fit.matrix.exp + error)
  results <- list(fit.matrix.exp=fit.matrix.exp, fit.matrix=fit.matrix)
  return(results)
}


#' Simulate data under additive model.
#'
#' @param n.muts Number of mutations.
#' @param addcoes Vector of additive effect for each mutation
#' @param sigma Noise to add to expectations
#' @param w.wt Fitness of the wildtype
#' @param geno.matrix Specifies which genotypes are simulated.
#' @return List:\cr
#' [[1]] \code{fit.matrix} is matrix with (simulated) observed data \cr
#' [[2]] \code{fit.matrix.exp} is same matrix without error
#' @details Function takes additive effects (\code{addcoes}) and generates the expected
#' fitness values for each genotype. This set of expected values is returned as \code{fit.matrix.exp}.
#' Then it adds normal error to them to produce the "observed" data (\code{fit.matrix}).
#' @examples
#' addcoes <- rep(0.3, 3)
#' geno.matrix <- generate.geno.matrix(3)
#' sim.add.data(3,addcoes,0.1,1,geno.matrix)
#' @export

sim.add.data <- function(n.muts, addcoes, sigma, w.wt, geno.matrix){

  n.genos <- length(geno.matrix[,1])
  geno.coes <- geno.matrix*addcoes
  #fit.matrix.exp <- as.matrix(apply(geno.coes, 1, function(x) w.wt + d.true*(1-prod(1-x))))
  fit.matrix.exp <- as.matrix(apply(geno.coes, 1, function(x) w.wt+sum(x)))
  error <- c(0, rnorm(n=n.genos-1, mean=0, sd=sigma))
  fit.matrix <- as.matrix(fit.matrix.exp + error)
  results <- list(fit.matrix.exp=fit.matrix.exp, fit.matrix=fit.matrix)
  return(results)
}



#' Simulate fitness data under multiplicative and additive models
#'
#' @param epi.model \code{mult/add} Epistasis model to simulate under.
#' @param mut.vals Vector of number of mutations to simulate
#' @param coe.vals Vector of stickbreaking coefficients to simulate
#' @param sig.vals Vector of sigma values to simulate
#' @param w.wt Fitness of the wildtype. Default 1.
#' @param n.reps.ea Number of replicates per parametric condition
#' @param print.status TRUE/FALSE. Should loop counters be printed.
#' @param outdir The path to write output files to (see details about file names).
#' @param wts Weight to give mutation on wildtype background vs other backgrounds. Default is c(2,1).
#'
#' @details Function contains a loop for combining each value of mut.vals, coe.vals and sig.vals,
#' generating data under the specified model and then fitting it. \cr
#'
#' Results are written to files; the name of the output files are formed by concatenating the outpath argument to
#' the epi.model argument. Separate files are generated for each
#' number of mutations (because the dimensionality of the output file changes with the number of mutations).
#' The output files are tab-delimited text files with one row per replicate.
#' The first 5 columns provide the parameter values and the rest of the columns give parameter estimates and
#' measures of fit.
#' @return Nothing. Instead results are written to output files and deposited in inst/extdata.
#' The files are named by appending the method
#' @examples
#' \dontrun{
#' sim.fit.mult.add.data.batch("mult",
#'   c(3,4,5),
#'   c(0.1, 0.3, 0.5),
#'   c(0.02, 0.05, 0.08),
#'   1,
#'   100,
#'   print.status=TRUE,
#'   outdir="~/Desktop",
#'   c(2,1))
#'   }
#' @export

sim.fit.mult.add.data.batch <- function(epi.model, mut.vals, coe.vals, sig.vals, w.wt=1, n.reps.ea=100, print.status=FALSE, outdir, wts){
    if (epi.model == "mult"){
      outfile <- "Mult_batch_out_test"
      outpath <- paste(outdir, outfile, sep="/")
    } else if (epi.model == "add"){
      outfile <- "Add_batch_out_test"
      outpath <- paste(outdir, outfile, sep="/")
    }
    for (mut.i in 1:length(mut.vals)){
    n.muts <- mut.vals[mut.i]
    geno.matrix <- generate.geno.matrix(n.muts)
    n.genos <- 2^n.muts
    if (print.status==TRUE) print(paste("n.muts =", n.muts))
    first.results <- TRUE

    for (coe.i in 1:length(coe.vals)){
      coes <- rep(coe.vals[coe.i], n.muts)
      if (print.status ==TRUE) print(paste("coe = ", coes[1]))

      for (sig.i in 1:length(sig.vals)){
        sigma <- sig.vals[sig.i]
        if (print.status == TRUE) print(paste("sigma = ", sigma))

        for (rep.i in 1:n.reps.ea){
          parm.row <- unlist(list(sim.model=epi.model, fit.model=epi.model, n.muts=n.muts, coes=coes[1], sigma=sigma, rep=rep.i))
          if (epi.model == "mult"){
            fit.matrix <- sim.mult.data(n.muts, selcoes=coes, sigma=sigma, w.wt=w.wt, geno.matrix)$fit.matrix
            fit.mult <- fit.mult.model(geno.matrix, fit.matrix, wts)
            fit.row <- unlist(list(round(unlist(fit.mult[1:4]),5), P=round(fit.mult$regression.results$P,5)))
          } else if (epi.model == "add"){
            fit.matrix <-  sim.add.data(n.muts, addcoes=coes, sigma=sigma, w.wt=w.wt, geno.matrix)$fit.matrix
            fit.add <- fit.add.model(geno.matrix, fit.matrix, wts)
            fit.row <- unlist(list(round(unlist(fit.add[1:4]),5), P=round(fit.add$regression.results$P,5)))
          }
          one.row <- as.data.frame(t(c(parm.row, fit.row)))
          output.file <- paste(outpath, "_", n.muts, "muts.txt", sep="")
          if (first.results == TRUE){
            write.table(x=one.row, file=output.file, append=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
            if (print.status ==TRUE)  print(paste("Writing results to: ", output.file), sep="")
            first.results <- FALSE
          } else{
            write.table(x=one.row, file=output.file, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE)
          }
        } #next rep.i
      } #next sig.i
    } # next coe.i
  } #next mut.i

}  #end simulate function
