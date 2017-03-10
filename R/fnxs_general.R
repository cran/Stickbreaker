
#' Internal simulation function to generate a matrix to weight the genotypes when estimating d and stickbreaking coefficients
#'
#' @param geno.matrix Genotype matrix generated in \code{\link{generate.geno.matrix}}
#' @param fit.matrix Fitness matrix generated in \code{\link{sim.stick.data}}
#' @param wts Vector of weights to be applied of form c(singletons, multiples). Default \code{wts=c(2,1)}.
#' @return weight.matrix
#' @details To calculate the likelihood of the data under the stickbreaking model for a given value of d,
#' we need to estimate the stickbreaking coefficients. The issue is whether all genotypes in the network provide equally good information
#' about the coefficients. The default assumption is that the wild type fitness
#' is know without error while all other genotypes have the same error structure. Coefficient estimates are based
#' on comparing pairs of genotypes: with and without the mutation. Therefore, estimates based on comparing wild type to
#' the single mutations (singletons) are expected to have half the variance of all other comparisons (i.e. multiples). This
#' function generates the weights matrix that reflects this. To change this assumption, change the \code{wts} parameter. For example,
#' if the wild type has the same error as all other genotypes, \code{wts = c(1,1)} would be appropriate.
#' @export

generate.geno.weight.matrix <- function(geno.matrix, fit.matrix, wts=c(2,1)){
  weight.matrix <- matrix(nrow=dim(geno.matrix)[1], ncol=dim(geno.matrix)[2])
  wt.row <- which(rowSums(geno.matrix)==0)
  mult.rows <- which(rowSums(geno.matrix)>1)
  single.rows <- which(rowSums(geno.matrix)==1)
  n.muts <- length(geno.matrix[1,])
  if (length(wts)==2){
    for (geno.i in 1:dim(geno.matrix)[1]){
      if (geno.i %in% single.rows){
        weight.matrix[geno.i,which(geno.matrix[geno.i,]==1)] <- wts[1]
      } else if (geno.i %in% mult.rows){
        weight.matrix[geno.i,which(geno.matrix[geno.i,]==1)] <- wts[2]
      }
    }
  } else{       # Using wts vector provided

    geno.sim.strings.h <- apply(geno.matrix, MARGIN=1, FUN=paste, collapse="")
    mut.i.genos <- apply(geno.matrix, 2, function(x) which(x==1))
    if (is.list(mut.i.genos)==FALSE){
      mut.i.genos <- as.list(as.data.frame(mut.i.genos))
    }
    for (mut.i in 1:n.muts){
      for (geno.i in 1:length(mut.i.genos[[mut.i]])){   # so geno.i is indexing mut.i.genos (not geno.matrix or fit.matrix)
        #geno.ii <- mut.i.genos[geno.i, mut.i]
        geno.ii <- mut.i.genos[[mut.i]][geno.i]
        geno <- geno.matrix[geno.ii,]
        geno.background <- geno
        geno.background[mut.i] <- 0
        geno.back.string <- paste(geno.background, collapse="")
        back.id <- which(geno.sim.strings.h==geno.back.string)
        var.of.diff <- wts[geno.ii] + wts[back.id]
        weight.of.geno <- 1/var.of.diff
        weight.matrix[geno.ii, mut.i] <- weight.of.geno
      } #next geno.i
    }  #next.mut.i

  }

  # Normalize columns (mutations) to sum to one
  for (mut.i in 1:dim(geno.matrix)[2]){
    weight.matrix[,mut.i] <- weight.matrix[,mut.i]/sum(weight.matrix[,mut.i], na.rm=T)
  }
  return(weight.matrix)
}

#' A dataset containing the fitness values for recombinant strains for Methylobacterium extorquens.
#' @references Chou, H., Chiu, H., Delaney, N., Segre, D., and Marx, C. 2011. Diminishing returns epistasis among beneficial mutations decelerates adaptation. Science 332, 1190-1192.
#' @usage data(Chou.data)
#' @format A list of 5 elements containing integers (1 or 0) that indicate the absence or presence of mutated alleles at 4 loci and a numerical element with relative fitness values.
#' @name Chou.data
NULL
#' A dataset containing the fitness values for recombinant poliovirus viruses.
#' @references Burns, C.C., Shaw, J., Campagnoli, R., Jorba, J., Vincent, A., Quay, J., and Kew, O. (2006). Modulation of Poliovirus Replicative Fitness in HeLa Cells by Deoptimization of Synonymous Codon Usage in the Capsid Region. J Virol 80, 3259-3272.
#' @usage data(burns.data)
#' @format A list of 5 elements containing integers (1 or 0) that indicate the absence or presence of mutated alleles at 4 loci and a numerical element with relative fitness values.
#' @name burns.data
NULL
#' A dataset containing the fitness values for recombinant Escherichia coli bacteria.
#' @references Khan, A. I., D. M. Dinh, D. Schneider, R. E. Lenski and T. F. Cooper, 2011 Negative epistasis between beneficial mutations in an evolving bacterial population. Science 332, 1193-1196.
#' @usage data(Khan.data)
#' @format A list of 6 elements containing integers (1 or 0) that indicate the absence or presence of mutated alleles at 5 loci and a numerical element with relative fitness values.
#' @name Khan.data
NULL
#' A dataset containing the fitness values for recombinant Escherichia coli bacteria.
#' @references Khan, A. I., D. M. Dinh, D. Schneider, R. E. Lenski and T. F. Cooper, 2011 Negative epistasis between beneficial mutations in an evolving bacterial population. Science 332, 1193-1196.
#' @usage data(Khan.data)
#' @format A list of 6 elements containing integers (1 or 0) that indicate the absence or presence of mutated alleles at 5 loci and a numerical element with relative fitness values.
#' @name caudle.data
NULL
