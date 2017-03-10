## ----setup, echo=FALSE, include=FALSE------------------------------------
suppressPackageStartupMessages({
  library(nnet)
  library(lattice)
  library(xtable)
  library(Stickbreaker)
})

## ---- echo=TRUE----------------------------------------------------------
  data <- Chou.data
  #data <- read.csv("your_data.csv")

## ---- echo=TRUE----------------------------------------------------------
  fit.results <- fit.models(data, d.range=c(0.1, 10), d.adj.max=1.1, wts=c(2,1))

## ---- echo=TRUE----------------------------------------------------------
  coes.prior <- c(0.05, 0.5)
  sig.prior <- c(0, 0.25)
  n.samps.per.mod <- 50   # set to smaller value when setting code up; larger value when doing real analysis
  d.range <- c(0.1, 10)
  d.adj.max <- 1.1
  wts <- c(2,1)
  min.R2 <- -1
  analysis.name <- "Test"
    
  posterior.results <- sim.data.calculate.posteriors(fit.results$fit.smry, data, analysis.name, coes.prior, sig.prior, d.range, d.adj.max, wts, n.samps.per.mod, min.R2, print.interval=25, tempdir())

## ------------------------------------------------------------------------
  data <- Khan.data
  #data <- read.csv("your_data.csv")

## ------------------------------------------------------------------------
  wts <- c(2,1)
  d.range <- c(0.1, 10)
  d.adj.max <- 1.1
  
  n.genos <- length(data[,1])
  n.muts <- length(data[1,])-1
  geno.matrix <- data[,seq(1, n.muts)]
  fit.matrix <- as.matrix(data[,(n.muts+1)])

## ------------------------------------------------------------------------
  d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
  d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
  d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, 
                                     d.hat.MLE, d.hat.RDB, d.range, d.adj.max)
  d.hat.MLE
  d.hat.RDB
  d.hat.seq

## ------------------------------------------------------------------------
  fit.stick <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, run.regression=TRUE)
  fit.stick

## ------------------------------------------------------------------------
  fit.mult <- fit.mult.model(geno.matrix, fit.matrix, wts=wts)
  #fit.mult

## ------------------------------------------------------------------------
  fit.add <- fit.add.model(geno.matrix, fit.matrix, wts=wts)
  #fit.add

## ------------------------------------------------------------------------
  outdir <- tempdir()
  plot.name <- "Ex_fit_vs_effect_plot.svg"
  file.path <- paste(outdir, plot.name, sep="/")
  svg(file=file.path, width=8, height=8)

  layout(mat=matrix(nrow=4, ncol=(1+n.muts), data=seq(1,4*(1+n.muts)), byrow=T), widths=c(0.25, rep(3,n.muts)), heights=c(0.25, rep(3,3))) -> l
  #layout.show(l)
  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  mod.names <- c("STICK", "MULT", "ADD")
  mut.names <- colnames(data)[1:n.muts]
  mod.cols <- c("deepskyblue","red", "yellow")
  for (i in 1:n.muts){
    plot(0,type='n',axes=FALSE,ann=FALSE, ylim=c(0,1), xlim=c(0,1))
    text(0.5, 0.5, labels=mut.names[i], font=2)
  }
  #--- STICK ---
  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE, ylim=c(0,1), xlim=c(0,1))
  text(0.5, 0.5, labels=mod.names[1], font=2, srt=90)
  par(mar=c(4,4,1,1))
  for (mut.i in 1:n.muts){
    plot(x=fit.stick$regression.results$fitness.of.backs[,mut.i], y=fit.stick$regression.results$effects.matrix[,mut.i], ylab="Effect", xlab="Back fitness", pch=21, bg=mod.cols[1])
    abline(fit.stick$regression.results$lm.intercepts[mut.i], fit.stick$regression.results$lm.slopes[mut.i])
  }
  #--- MULT ---
  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE, ylim=c(0,1), xlim=c(0,1))
  text(0.5, 0.5, labels=mod.names[2], font=2, srt=90)
  par(mar=c(4,4,1,1))
  for (mut.i in 1:n.muts){
    plot(x=fit.mult$regression.results$fitness.of.backs[,mut.i], y=fit.mult$regression.results$effects.matrix[,mut.i], ylab="Effect", xlab="Back fitness", pch=21, bg=mod.cols[2])
    abline(fit.mult$regression.results$lm.intercepts[mut.i], fit.mult$regression.results$lm.slopes[mut.i])
  }
  #--- ADD ---
  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE, ylim=c(0,1), xlim=c(0,1))
  text(0.5, 0.5, labels=mod.names[2], font=2, srt=90)
  par(mar=c(4,4,1,1))
  for (mut.i in 1:n.muts){
    plot(x=fit.mult$regression.results$fitness.of.backs[,mut.i], y=fit.mult$regression.results$effects.matrix[,mut.i], ylab="Effect", xlab="Back fitness", pch=21, bg=mod.cols[2])
    abline(fit.mult$regression.results$lm.intercepts[mut.i], fit.mult$regression.results$lm.slopes[mut.i])
  }
  dev.off()
  

## ---- out.width = "800px"------------------------------------------------
    knitr::include_graphics(file.path)

## ------------------------------------------------------------------------
  outdir <- tempdir()
  plot.name <- "Ex_pred_obs_fit.svg"
  file.path <- paste(outdir, plot.name, sep="/")
  svg(file=file.path, width=8, height=6)

  mod.cols <- c("deepskyblue","red", "yellow")
  lims <- c(1,1.4)
  par(mar=c(5,4,4,1))
  plot(x=fit.matrix, y=fit.stick$pred.matrix$pred, ylim=lims, xlim=lims, ylab="Model Predicted Fitness", xlab="Observed Fitness", pch=21, bg=mod.cols[1], cex=1.1)
  abline(0,1, lty="dashed")
  points(x=fit.matrix, y=fit.mult$pred.matrix$pred, pch=21, bg=mod.cols[2], cex=0.9)  
  points(x=fit.matrix, y=fit.add$pred.matrix$pred, pch=21, bg=mod.cols[3], cex=0.8) 
  text(x=fit.matrix, y=fit.stick$pred.matrix$pred, labels=fit.stick$pred.matrix$string, srt=90, cex=0.7, pos=3, off=2)
  legend("topleft", legend=c("Stick", "Mult", "Add"), pch=21, pt.bg=mod.cols, bty="n")
  dev.off()

## ---- out.width = "800px"------------------------------------------------
    knitr::include_graphics(file.path)

## ------------------------------------------------------------------------
  n.muts <- 3
  coe <- 0.1
  coes <- rep(coe, n.muts)
  sigma <- 0.1
  w.wt <- 1         #fitness of wild type
  d.true <- 1       # distance to fitness boundary
  n.genos <- 2^n.muts     # number of genotypes in full network
  d.range <- c(0.1, 10)
  d.adj.max <- 1.1

  geno.matrix <- generate.geno.matrix(n=n.muts)
  stick.sim.data <- sim.stick.data(n.muts, coes, sigma, d.true, w.wt, geno.matrix)
  fit.matrix <- stick.sim.data$fit.matrix
  #print(fit.matrix)

## ------------------------------------------------------------------------
  n.muts <- 3
  selcoe <- 0.3
  selcoes <- rep(selcoe, n.muts)
  sigma <- 0.1
  w.wt <- 1         #fitness of wild type
  n.genos <- 2^n.muts     # number of genotypes in full network

  geno.matrix <- generate.geno.matrix(n=n.muts)
  mult.sim.data <- sim.mult.data(n.muts, selcoes, sigma, w.wt, geno.matrix)
  fit.matrix <- mult.sim.data$fit.matrix
  #print(fit.matrix)

## ------------------------------------------------------------------------
  n.muts <- 3
  addcoe <- 0.3
  addcoes <- rep(addcoe, n.muts)
  sigma <- 0.1
  w.wt <- 1         #fitness of wild type
  n.genos <- 2^n.muts     # number of genotypes in full network

  geno.matrix <- generate.geno.matrix(n=n.muts)
  add.sim.data <- sim.add.data(n.muts, addcoes, sigma, w.wt, geno.matrix)
  fit.matrix <- add.sim.data$fit.matrix
  #print(fit.matrix)

## ------------------------------------------------------------------------
  sim.batch.data <- TRUE
  
  mut.vals <- c(3,4,5)
  coe.vals <- c(0.1, 0.3, 0.5)
  sig.vals <- c(0.02, 0.05, 0.08)
  w.wt <- 1
  d.true <- 1
  wts <- c(2,1)
  n.reps.ea <- 10  # set to a small value just to make code run fast 
  d.range <- c(0.1, 10)
  fit.methods <- "seq"
  d.max.adj <- 1.0
  run.regression <- FALSE  
    
  if (sim.batch.data == TRUE){
    outdir <- tempdir()
    sim.fit.stick.data.batch(mut.vals=mut.vals, coe.vals=coe.vals, sig.vals=sig.vals, d.true=d.true, d.range=d.range, w.wt=w.wt, n.reps.ea=n.reps.ea, print.status=FALSE, fit.methods=fit.methods, outdir=outdir, wts, d.max.adj=d.max.adj, run.regression, RDB.method="pos")
  }

## ------------------------------------------------------------------------
  sim.batch.data <- TRUE
  
  mut.vals <- c(3,4,5)
  coe.vals <- c(0.1, 0.3, 0.5)
  sig.vals <- c(0.02, 0.05, 0.08)
  w.wt <- 1
  n.reps.ea <- 10
  
  #--- mulitiplicative ---
  if (sim.batch.data == TRUE){
    outdir <- tempdir()
    sim.fit.mult.add.data.batch(epi.model="mult", mut.vals=mut.vals, coe.vals=coe.vals, sig.vals=sig.vals, w.wt=w.wt, n.reps.ea=n.reps.ea, print.status=FALSE, outdir=outdir, wts)
  }

## ------------------------------------------------------------------------
  sim.batch.data <- TRUE
  
  mut.vals <- c(3,4,5)
  coe.vals <- c(0.1, 0.3, 0.5)
  sig.vals <- c(0.02, 0.05, 0.08)
  w.wt <- 1
  n.reps.ea <- 10
  wts <- c(2,1)
  
  if (sim.batch.data == TRUE){
    outdir <- tempdir()
    sim.fit.mult.add.data.batch(epi.model="add", mut.vals=mut.vals, coe.vals=coe.vals, sig.vals=sig.vals, w.wt=w.wt, n.reps.ea=n.reps.ea, print.status=FALSE, outdir=outdir, wts)
  }

## ------------------------------------------------------------------------

  sim.training.data.from.priors <- TRUE

  if (sim.training.data.from.priors == TRUE){
    coes.prior <- c(0.05, 0.5)
    sig.prior <- c(0, 0.25)
    n.samps.per.mod <- 50
    d.true <- 1
    d.range <- c(0.1, 10)
    d.adj.max <- 1.1
    w.wt <- 1
    wts <- c(2,1)
    muts.to.sim <- c(4)   # vector; to simulate, for example, sets of 3, 4 and 5 mutations at once, use c(3,4,5)
    print.interval <- NA  # set to NA to block update printing
    
    for (mut.i in 1:length(muts.to.sim)){
      n.muts <- muts.to.sim[mut.i]
      print(paste("n.muts=", n.muts))
      outdir <- tempdir()
      file.name <- paste("Training_simulated_priors_fit_data_", n.muts, "muts.txt", sep="")
      outpath <- paste(outdir, file.name, sep="/")
      sim.data.from.priors.for.mod.selection(n.muts=n.muts, coes.prior=coes.prior, sigs.prior=sig.prior, mods.to.sim=c("stick", "mult", "add"), d.true=d.true, d.range=d.range, d.adj.max=d.adj.max, w.wt=w.wt, wts=wts, outpath=outpath, n.samps.per.mod=n.samps.per.mod, coe.sim.model="identical", coe.dist.par=NA, print.interval=print.interval)
    } #next mut.i
  }

## ---- eval=TRUE, results="hide"------------------------------------------
  
  model.file <- "nnet_mod_R2_P"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  data.file <- "Training_simulated_priors_fit_data"

  fit.multinomial <- TRUE
  min.R2 <- -1
  
  if (fit.multinomial == TRUE){
  
    mod.formula <- as.formula(model ~ R2.stick + R2.mult + R2.add + P.stick + P.mult + P.add)
    muts.to.sim <- c(4)
    fit.nmet.models <- vector("list", length(muts.to.sim))
    for (mut.i in 1:length(muts.to.sim)){
      n.muts <- muts.to.sim[mut.i]
      indir <- tempdir()
      file.name <- paste(data.file, "_", n.muts, "muts.txt", sep="")
      inpath <- paste(indir, file.name, sep="/")
      if (file.exists(inpath)){
        data <- read.table(file=inpath, header=TRUE)
        data$R2.stick[which(data$R2.stick < min.R2)] <- min.R2
        data$R2.mult[which(data$R2.mult < min.R2)] <- min.R2
        data$R2.add[which(data$R2.add < min.R2)] <- min.R2
        fit.nmet.models[[mut.i]] <- fit.nnet.multinomial.regression(data, mod.formula)
        m <- fit.nnet.multinomial.regression(data, mod.formula)
        modelout <- paste(modpath, "_", n.muts, "muts.rda", sep="")
        saveRDS(m, file=modelout)
      }
    }  
  }

## ------------------------------------------------------------------------
  n.muts <- 5
  coes <- 0.2
  sigma <- 0.03
  geno.matrix <- generate.geno.matrix(5)
  geno.matrix <- geno.matrix[which(rowSums(geno.matrix) %in% c(0,1,4,5)),]
  fit.matrix <- sim.stick.data(n.muts, coes=coes, sigma=sigma, d.true=1, w.wt=1, geno.matrix)

## ---- eval=TRUE----------------------------------------------------------
  outdir <- tempdir()
  file.name <- paste("Training_simulated_priors_fit_data_first_last", n.muts, "muts.txt", sep="")
  simdata.outpath <- paste(outdir, file.name, sep="/")
  coes.prior <- c(0.05, 0.5)
  sig.prior <- c(0, 0.25)
  n.samps.per.mod <- 10   # Increase to a large number when analyzing real data
  d.true <- 1
  d.range <- c(0.1, 10)
  d.adj.max <- 1.1
  w.wt <- 1
  wts <- c(2,1)
  n.samps.per.mod <- 50

  sim.partial.data.from.priors.for.mod.selection(geno.matrix,
                                                      coes.prior=coes.prior,
                                                      sigs.prior=sig.prior,
                                                      mods.to.sim=c("stick", "mult", "add"),
                                                      d.true=1,
                                                      d.range=d.range,
                                                      d.adj.max=d.adj.max,
                                                      w.wt=1,
                                                      wts=wts,
                                                      outpath=simdata.outpath,
                                                      n.samps.per.mod=n.samps.per.mod,
                                                      coe.sim.model="identical",
                                                      coe.dist.par=NA,
                                                      print.interval=25)


## ---- eval=TRUE, results="hide"------------------------------------------
  
  model.file <- "nnet_mod_R2_P_5muts_first_last"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")

  fit.multinomial <- TRUE
  min.R2 <- -1
  
  if (fit.multinomial == TRUE){
    mod.formula <- as.formula(model ~ R2.stick + R2.mult + R2.add + P.stick + P.mult + P.add)
    indir <- tempdir()
    file.name <- paste("Training_simulated_priors_fit_data_first_last", n.muts, "muts.txt", sep="")
    inpath <- paste(indir, file.name, sep="/")
      if (file.exists(inpath)){
        data <- read.table(file=inpath, header=TRUE)
        data$R2.stick[which(data$R2.stick < min.R2)] <- min.R2
        data$R2.mult[which(data$R2.mult < min.R2)] <- min.R2
        data$R2.add[which(data$R2.add < min.R2)] <- min.R2
        m <- fit.nnet.multinomial.regression(data, mod.formula)
        modelout <- paste(modpath, "_", n.muts, "muts.rda", sep="")
        saveRDS(m, file=modelout)
      }
    }  

## ------------------------------------------------------------------------

  sim.testing.data.from.vectors <- TRUE
  testdata.file <- "Test_simulated_fit_data"
  

  if (sim.testing.data.from.vectors == TRUE){
    mods.to.sim <- c("stick", "mult", "add")   # only the strings "stick", "mult" and "add" are allowed
    muts.to.sim <- c(4)   # c(3,4,5)
    coes.to.sim <- c(0.1, 0.3)   # c(0.1, 0.2, 0.3, 0.4, 0.5)
    sigs.to.sim <- c(0.03, 0.06)
    d.true <- 1
    d.range <- c(0.1, 10)
    w.wt <- 1
    wts <- c(2,1)
    n.reps.ea <- 100
    for (mut.i in 1:length(muts.to.sim)){
      n.muts <- muts.to.sim[mut.i]
      print(paste("n.muts=", n.muts))
      outdir <- tempdir()
      file.name <- paste(testdata.file, "_", n.muts, "muts.txt", sep="")
      outpath <- paste(outdir, file.name, sep="/")
      sim.data.for.mod.selection(n.muts=n.muts, coes.to.sim, sigs.to.sim, mods.to.sim, d.true, d.range, w.wt,wts,  outpath, n.reps.ea=n.reps.ea)
    }
  }

## ------------------------------------------------------------------------
  fit.test.data <- TRUE

  model.file <- "nnet_mod_R2_P"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  testdata.file <- "Test_simulated_fit_data"

  muts.to.lyze <- c(4)   # The mutations to analyze
  
  
  if (fit.test.data == TRUE){
  
    for (mut.i in 1:length(muts.to.lyze)){
      n.muts <- muts.to.lyze[mut.i]
      print(paste("muts =", n.muts))
      modelin <- paste(modpath, "_", n.muts, "muts.rda", sep="")
      mod <- readRDS(modelin)
      
      datadir <- tempdir()
      file.name <- paste(datadir, "/", testdata.file, "_", n.muts, "muts.txt", sep="")
      test.data <- read.table(file.name, header=TRUE)
    
      posteriors <- as.data.frame(matrix(nrow=length(test.data[,1]), ncol=3))
      colnames(posteriors)=c("add", "mult", "stick")
    
      for (j in 1:length(test.data[,1])){
        posteriors[j,] <- predict(mod, newdata=test.data[j,], type="probs")
      }
    
      test.data <- cbind(test.data, posteriors)  
      outdata.file <- file.name <- paste(datadir, "/", testdata.file, "_", n.muts, "muts_wPosts.csv", sep="")
      write.csv(x=test.data, file=file.name, row.names=FALSE)
    }  #next mut.i
  }  # end if fit.test.data==TRUE

## ------------------------------------------------------------------------
 testdata.file <- "Test_simulated_fit_data"
  muts.to.lyze <- c(4)
  rej.cut <- 0.05       # consider model rejected when it has posterior < this
  posterior.list <- vector("list", length(muts.to.lyze))
  names(posterior.list) <- paste("muts.",  muts.to.lyze, sep="")
  
  for (mut.i in 1:length(muts.to.lyze)){
     n.muts <- muts.to.lyze[mut.i]
     datadir <- tempdir()
     indata.file <-  paste(datadir, "/", testdata.file, "_", n.muts=4, "muts_wPosts.csv", sep="")
     data <- read.csv(file=indata.file, stringsAsFactors = FALSE)
  
     posterior.list[[mut.i]] <- summarize.posteriors.on.simulated.dataset(data, rej.cut, n.muts)
  }
posterior.list[[1]]   # print to illustrate output

## ------------------------------------------------------------------------
	data(caudle.data)
	fit.col <- which(colnames(caudle.data)=="Fitness")
	n.singles <- fit.col-1
	n.muts <- n.singles

	fitness.shift <- -caudle.data$Fitness[1] + 1      # add this much to each fitness. Need wt to have positive fitness
	caudle.data$Fitness <- caudle.data$Fitness + fitness.shift
	geno.matrix <- geno.matrix.Caud <- caudle.data[,1:(fit.col-1)]
	fit.matrix <- fit.matrix.Caud <- t(t(caudle.data[,fit.col:fit.col]))
	wts <- c(2,1)
	d.range <- c(max(fit.matrix)-fit.matrix[1], 2*max(fit.matrix))
  d.adj.max <- 1.1
  d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range, wts)
  d.vect <- seq(d.range[1], d.range[2], length.out=100)
  
  # Visualize the likelihood of d
  like.profile <- calc.stick.logLn(geno.matrix, fit.matrix, d.vect, wts)
  plot(d.vect, like.profile, type="l", xlim=c(15,35), ylab="Log-likelihood", xlab="d")
  abline(v=d.hat.MLE, lty="dashed")
  
  d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix, no.est=NA)$d.hat.RDB  # doesn't work for double data
  d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range)

  fit.stick <- fit.stick.Caud <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, wts=wts, run.regression=TRUE)
  fit.mult <- fit.mult.Caud <- fit.mult.model(geno.matrix, fit.matrix, wts=wts)
  fit.add <- fit.add.Caud <- fit.add.model(geno.matrix, fit.matrix, wts=wts)
  fit.smry <- fit.smry.Caud <- summarize.fits.for.posterior.calc(fit.stick, fit.mult, fit.add)
  
  nmuts <- apply(geno.matrix, 1, sum)
  col.v <- c("deepskyblue", "red", "yellow")

  plot(x=fit.matrix, y=fit.stick$pred.matrix$pred, ylim=c(0,20), xlim=c(0,20), ylab="Model predicted fitness", xlab="Rescaled observed fitness", pch=21, bg=col.v[1], cex=1.1)
  abline(0,1, lty="dashed")
  points(x=fit.matrix, y=fit.mult$pred.matrix$pred, pch=21, bg=col.v[2], cex=0.9)  
  points(x=fit.matrix, y=fit.add$pred.matrix$pred, pch=21, bg=col.v[3], cex=0.8)  
	legend("topleft", pch=21, pt.bg=col.v, bty="n", legend=c("Stick", "Mult", "Add"))

## ------------------------------------------------------------------------

  sim.training.data.from.priors <- TRUE

  if (sim.training.data.from.priors == TRUE){
    mods.to.sim <- c("stick", "mult", "add")
    coes.prior <- c(0.05, 0.5)
    sigs.prior <- c(0, 0.25)
    n.samps.per.mod <- 50
    
    d.true <- 1
    d.range <- c(0.1, 10)
    d.adj.max <- 1.1
    w.wt <- 1
    wts <- c(2,1)
    print.interval <- 25
    
    outdir <- tempdir()
    file.name <- paste("Training_simulated_priors_fit_caudle_data.txt", sep="")
    outpath <- paste(outdir, file.name, sep="/")
    
    sim.partial.data.from.priors.for.mod.selection(geno.matrix, coes.prior=coes.prior, sigs.prior=sigs.prior, mods.to.sim=mods.to.sim, d.true=d.true, d.range=d.range, d.adj.max=d.adj.max, w.wt=w.wt, wts=wts, outpath=outpath, n.samps.per.mod=n.samps.per.mod, coe.sim.model="identical", coe.dist.par=NA, print.interval=25)
  }

## ---- eval=TRUE, results="hide"------------------------------------------
  
  model.file <- "nnet_mod_R2_P_caudle"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  data.file <- "Training_simulated_priors_fit_caudle_data.txt"
  min.R2 <- -1

  fit.multinomial <- TRUE
  
  if (fit.multinomial == TRUE){
  
    mod.formula <- as.formula(model ~ R2.stick + R2.mult + R2.add + P.stick + P.mult + P.add)
    indir <- tempdir()
    inpath <- paste(indir, data.file, sep="/")
    if (file.exists(inpath)){
      data <- read.table(file=inpath, header=TRUE)
      data$R2.stick[which(data$R2.stick<min.R2)] <- min.R2
      data$R2.mult[which(data$R2.mult<min.R2)] <- min.R2
      data$R2.add[which(data$R2.add<min.R2)] <- min.R2
      fit.nmet.model <- fit.nnet.multinomial.regression(data, mod.formula)
      m <- fit.nnet.multinomial.regression(data, mod.formula)
      modelout <- paste(modpath, ".rda", sep="")
      saveRDS(m, file=modelout)
    }
  }

## ---- eval=TRUE----------------------------------------------------------
  model.file <- "nnet_mod_R2_P_caudle"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  modelin <- paste(modpath, ".rda", sep="")
  model <- readRDS(modelin)
  post.data <- post.data.Caud <- calculate.posteriors.for.datasets(fit.smry.Caud,  model)
  smry.w.probs.Caud <- as.list(post.data)
  print(post.data.Caud)


## ------------------------------------------------------------------------
  #datadir <- system.file("data", package="Stickbreaker")
  #load(paste(datadir, "burns.data.rda", sep="/"))
  data <- burns.data
  wts <- c(2,1)
  n.genos <- length(data[,1])
  n.muts <- length(data[1,])-1
  mut.i <- n.muts - 2
  n.mut.i <- mut.i
  geno.matrix <- geno.matrix.Burns <- data[,seq(1, n.muts)]
  fit.matrix <- fit.matrix.Burns <-  as.matrix(data[,(n.muts+1)])
  fit.matrix <- log(fit.matrix)
  muts.by.geno <- apply(geno.matrix, 1, sum)
  
  d.range <- c(min(fit.matrix), max(fit.matrix)*1.5)
  d.adj.max <- 1.1
  d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix, d.range=d.range)
  d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix)$d.hat.RDB
  d.hat.seq <- estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range)
  fit.stick <- fit.stick.Burns <- fit.stick.model.given.d(geno.matrix, fit.matrix, d.hat.seq, run.regression=TRUE)
  fit.mult <- fit.mult.Burns <- fit.mult.model(geno.matrix, fit.matrix, wts=wts)
  fit.add <- fit.add.Burns <- fit.add.model(geno.matrix, fit.matrix, wts=wts)
  fit.smry <- fit.smry.Burns <- summarize.fits.for.posterior.calc(fit.stick, fit.mult, fit.add)
  
  lims=c(min(fit.matrix)*0.95, max(fit.matrix)*1.05)
  preds <- as.data.frame(cbind(obs=fit.stick$pred.matrix$fit,stick=fit.stick$pred.matrix$pred, mult=fit.mult$pred.matrix$pred, add=fit.add$pred.matrix$pred))
  rownames(preds) <- fit.stick$pred.matrix$string
  max.fits <- apply(preds, 1, function(x) max(x[2:4]))
  min.fits <- apply(preds, 1, function(x) min(x[2:4]))
  
  cex.v <- c(1,1,1)
  lab.cex <- 0.75
  leg.location <- "topleft"
  leg.cex <- 1.2
  col.v <- c("deepskyblue", "red", "yellow")
  plot(y=preds$stick, x=preds$obs, pch=21, bg=col.v[1], xlim=lims, ylim=lims, cex=cex.v[1], xlab="Rescaled observed fitness", ylab="Model predicted fitness")
  abline(0,1)
  points(y=preds$mult, x=preds$obs, pch=21, bg=col.v[2], cex=cex.v[2])
  points(y=preds$add, x=preds$obs, pch=21, bg=col.v[3], cex=cex.v[3])
  text(y=min.fits, x=preds$obs, pos=1, off=1, srt=90, labels=rownames(preds), cex=lab.cex)
  legend("topleft", bty="n", pch=rep(21, 3), pt.bg=col.v, legend=c("Stick", "Mult", "Add"), cex=leg.cex)

## ---- message=FALSE, warning=FALSE---------------------------------------

  sim.training.data.from.priors <- TRUE

  if (sim.training.data.from.priors == TRUE){
    mods.to.sim <- c("stick", "mult", "add")
    coes.prior <- c(-0.05, -1)
    sigs.prior <- c(0, 0.25)
    n.samps.per.mod <- 50
    
    d.true <- 1
    d.range <- c(1, 10)
    d.adj.max <- 1.1
    w.wt <- fit.matrix[1,1]
    wts <- c(2,1)
    print.interval <- 25
    
    outdir <- tempdir()
    file.name <- paste("Training_simulated_priors_fit_burns2_data.txt", sep="")
    outpath <- paste(outdir, file.name, sep="/")
    
    sim.partial.data.from.priors.for.mod.selection(geno.matrix, coes.prior=coes.prior, sigs.prior=sigs.prior, mods.to.sim=mods.to.sim, d.true=d.true, d.range=d.range, d.adj.max, w.wt=w.wt, wts=wts, outpath=outpath, n.samps.per.mod=n.samps.per.mod, coe.sim.model="identical", coe.dist.par=NA, print.interval=print.interval)

  }

## ---- eval=TRUE, results="hide"------------------------------------------
  
  model.file <- "nnet_mod_R2_P_burns2"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  data.file <- "Training_simulated_priors_fit_burns2_data.txt"

  fit.multinomial <- TRUE
  min.R2 <- -1
  
  if (fit.multinomial == TRUE){
  
    mod.formula <- as.formula(model ~ R2.stick + R2.mult + R2.add + P.stick + P.mult + P.add)
    indir <- tempdir()
    inpath <- paste(indir, data.file, sep="/")
    if (file.exists(inpath)){
      data <- read.table(file=inpath, header=TRUE)
      data$R2.stick[which(data$R2.stick<min.R2)] <- min.R2
      data$R2.mult[which(data$R2.mult<min.R2)] <- min.R2
      data$R2.add[which(data$R2.add<min.R2)] <- min.R2
      fit.nmet.model <- fit.nnet.multinomial.regression(data, mod.formula)
      m <- fit.nnet.multinomial.regression(data, mod.formula)
      modelout <- paste(modpath, ".rda", sep="")
      saveRDS(m, file=modelout)
    }
  }

## ------------------------------------------------------------------------
  model.file <- "nnet_mod_R2_P_burns2"
  moddir <- tempdir()
  modpath <- paste(moddir, model.file, sep="/")
  modelin <- paste(modpath, ".rda", sep="")
  mod <- readRDS(modelin)
  post.probs <- post.probs.Burns <- calculate.posteriors.for.datasets(fit.smry.Burns, mod=mod)
  print(post.probs.Burns)

