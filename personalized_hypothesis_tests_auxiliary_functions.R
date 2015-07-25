#######################################
## Auxiliary functions for generating
## the paper figures
#######################################


#######################################
## Auxiliary functions for Figure 3
#######################################

PlotBeforeAfterMedication <- function(healthCodeId, 
                                      featName, 
                                      dat, 
                                      ylab = "",
                                      xlab = "date",
                                      ylim = NULL,
                                      cex.dot = 1.5,
                                      legendPos = "topright",
                                      main = "",
                                      cex.lab = 1,
                                      cex.axis = 1,
                                      cex.main = 1,
                                      lcol = "black") {
  idx <- which(dat$healthCode == healthCodeId)
  sdat <- dat[idx, ]
  sdat <- sdat[order(sdat$createdOn),]
  ibefore <- which(sdat$momentInDayFormat.json.choiceAnswers == "Immediately before Parkinson medication")
  iafter <- which(sdat$momentInDayFormat.json.choiceAnswers == "Just after Parkinson medication (at your best)")
  feat <- sdat[, featName]
  plot(sdat$createdOn, feat, type = "l", main = main, xlab = xlab, ylab = ylab, 
       ylim = ylim, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, col = lcol)
  legend(legendPos, legend = c("before", "after", "other / missing"), 
         text.col = c("red", "blue", "grey"),
         bty = "n")
  points(sdat$createdOn, feat, col = "grey", pch = 20, cex = cex.dot)
  points(sdat$createdOn[ibefore], feat[ibefore], col = "red", pch = 20, cex = cex.dot)
  points(sdat$createdOn[iafter], feat[iafter], col = "blue", pch = 20, cex = cex.dot) 
}

ShuffleBeforeAfterLabels <- function(dat, seed) {
  ibefore <- which(dat$momentInDayFormat.json.choiceAnswers == "Immediately before Parkinson medication")
  iafter <- which(dat$momentInDayFormat.json.choiceAnswers == "Just after Parkinson medication (at your best)")
  nb <- length(ibefore)
  na <- length(iafter)
  ii <- c(ibefore, iafter)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  iib <- sample(ii, nb, replace = FALSE)
  iia <- setdiff(ii, iib)
  sdat <- dat
  sdat[iib, "momentInDayFormat.json.choiceAnswers"] <- "Immediately before Parkinson medication"
  sdat[iia, "momentInDayFormat.json.choiceAnswers"] <- "Just after Parkinson medication (at your best)"
  sdat
}

TTestPerm <- function(B, x1, x2) {
  TTestStat <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    s2.1 <- var(x1)
    s2.2 <- var(x2)
    xbar.1 <- mean(x1)
    xbar.2 <- mean(x2)
    s2 <- ((n1 - 1) * s2.1 + (n2 - 1) * s2.2)/(n1 + n2 - 2)
    (xbar.1 - xbar.2)/sqrt(s2 * (n1^(-1) + n2^(-1)))
  }
  n1 <- length(x1)
  n2 <- length(x2)
  x <- c(x1, x2)
  n <- n1 + n2
  i1 <- seq(n1)
  i2 <- seq(n1 + 1, n1 + n2, by = 1)
  stat <- rep(NA, B) 
  for (i in seq(B)) {
    xperm <- x[sample(n, replace = FALSE)]
    stat[i] <- TTestStat(xperm[i1], xperm[i2])
  }
  obs <- TTestStat(x1, x2)
  aobs <- abs(obs)
  pval <- (sum(stat >= aobs) + sum(stat <= -aobs))/B
  list(pval = pval, stat = stat)
}

WilcoxTestPerm <- function(B, x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  x <- c(x1, x2)
  n <- n1 + n2
  i1 <- seq(n1)
  i2 <- seq(n1 + 1, n1 + n2, by = 1)
  stat <- rep(NA, B) 
  for (i in seq(B)) {
    xperm <- x[sample(n, replace = FALSE)]
    stat[i] <- wilcox.test(xperm[i1], xperm[i2])$statistic
  }
  obs <- wilcox.test(x1, x2)$statistic
  aobs <- abs(obs)
  pval <- (sum(stat >= aobs) + sum(stat <= -aobs))/B
  list(pval = pval, stat = stat)
}

GetNullModelAnalyticalPvalues <- function(nSim, dat) {
  GenerateNullData <- function(dat, ib, ia) {
    dat$momentInDayFormat.json.choiceAnswers[ib] <- "Immediately before Parkinson medication"
    dat$momentInDayFormat.json.choiceAnswers[ia] <- "Just after Parkinson medication (at your best)"
    dat
  }  
  pvalsT <- rep(NA, nSim)
  pvalsW <- rep(NA, nSim)
  n <- nrow(dat)
  for (i in seq(nSim)) {
    cat(i, "\n")
    nb <- sample(seq(10, n-10, by = 1), 1)
    ib <- sample(seq(n), nb, replace = FALSE)
    ia <- setdiff(seq(n), ib)
    nulldat <- GenerateNullData(dat, ib, ia)
    xa <- nulldat[which(nulldat$momentInDayFormat.json.choiceAnswers == "Just after Parkinson medication (at your best)"), "numberTaps"]
    xb <- nulldat[which(nulldat$momentInDayFormat.json.choiceAnswers == "Immediately before Parkinson medication"), "numberTaps"]
    pvalsT[i] <- t.test(xa, xb, var.equal = TRUE)$p.value
    pvalsW[i] <- wilcox.test(xa, xb)$p.value
  }
  list(pvalsT = pvalsT, pvalsW = pvalsW)
}


#######################################
## auxiliary functions for Figure 4
#######################################

CheckNullUIPvalDistribution <- function(dat, 
                                        participantId, 
                                        nperm = 1000, 
                                        breakFeatCor = TRUE) {
  pdat <- GetParticipantBeforeAfterData(dat, participantId)[, 8:32]
  n <- nrow(pdat)
  featnms <- colnames(pdat)[-1]
  nfeat <- length(featnms)
  pttest <- rep(NA, nperm)
  pwilcoxson <- rep(NA, nperm)
  pvalsTtestM <- matrix(NA, nperm, nfeat)
  pvalsWilcoxM <- matrix(NA, nperm, nfeat)
  colnames(pvalsTtestM) <- featnms
  colnames(pvalsWilcoxM) <- featnms
  respName <- "momentInDayFormat.json.choiceAnswers"
  for (i in seq(nperm)) {
    cat(i, "\n")
    pdat2 <- pdat
    if (breakFeatCor) {
      for (j in seq(ncol(pdat))) {
        pdat2[, j] <- pdat2[sample(n), j]
      }
    }
    else {
      pdat2[, respName] <- pdat2[sample(n), respName] 
    }
    aux <- UITests(pdat2, mtmethod = "BH", alternative = "two.sided", respName)
    pttest[i] <- aux$pvalTtest
    pwilcoxson[i] <- aux$pvalWilcox
    pvalsTtestM[i,] <- aux$pvalsTtest
    pvalsWilcoxM[i,] <- aux$pvalsWilcox
  }
  list(pttest = pttest, 
       pwilcoxson = pwilcoxson,
       pvalsTtestM = pvalsTtestM, 
       pvalsWilcoxM = pvalsWilcoxM)
}


#######################################
## auxiliary functions for Figure 5
#######################################

RandomForestPermNull <- function(nperm, 
                                 dat,
                                 nSplits,
                                 splitSeed = NULL,
                                 respName, 
                                 refClass = "before") {
  respColumn <- which(colnames(dat) == respName)
  paucs <- rep(NA, nperm)
  if (!is.null(splitSeed)) {
    set.seed(splitSeed)
  }
  splitDat <- TrainTestRandomSplit(dat, nSplits)
  trainDat <- RemoveNAs(splitDat$trainDat)
  testDat <- RemoveNAs(splitDat$testDat)
  myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
  fit <- randomForest(myFormula, data = trainDat)
  pred <- predict(fit, testDat[, -respColumn], type = "prob")[, refClass]
  ties <- GetTieStats(pred)
  ytest <- testDat[, respColumn]
  nb <- sum(ytest == "before")
  na <- sum(ytest == "after")
  n <- length(ytest)
  for (i in seq(nperm)) {
    cat(i, "\n")
    predobj <- prediction(pred, ytest[sample(n, replace = FALSE)])
    paucs[i] <- performance(predobj, "auc")@y.values[[1]]
  }
  list(paucs = paucs, nb = nb, na = na, aux = ties$aux) 
}


GetNullModelAnalyticalPvalues2 <- function(nSim, dat) {
  GenerateNullData <- function(dat, ib, ia) {
    dat$momentInDayFormat.json.choiceAnswers[ib] <- "before"
    dat$momentInDayFormat.json.choiceAnswers[ia] <- "after"
    dat
  }  
  pvalsRF <- rep(NA, nSim)
  n <- nrow(dat)
  for (i in seq(nSim)) {
    cat(i, "\n")
    nb <- sample(seq(10, n-10, by = 1), 1)
    ib <- sample(seq(n), nb, replace = FALSE)
    ia <- setdiff(seq(n), ib)
    nulldat <- GenerateNullData(dat, ib, ia)
    ## we need a try here because due to the sometimes
    ## small nb (could be 10) when we do the random split
    ## into train/test sets we can get a unique label type
    ## in the test set
    aux <- try(RandomForestTest(nRuns = 1, nulldat, nSplits = 2, splitSeeds = 10000 + i, 
                                respName = "momentInDayFormat.json.choiceAnswers", 
                                refClass = "before")$pval, silent = TRUE)
    if (!inherits(aux, "try-error")) {
      pvalsRF[i] <- aux
    }
  }
  pvalsRF
}


#######################################
## auxiliary functions for Figure 7
#######################################

GetEmpiricalPower <- function(x, alphaGrid = seq(0.001, 1, by = 0.001)) {
  EmpiricalPower <- function(x, alphaGrid) {
    ntests <- length(x)
    nalpha <- length(alphaGrid)
    out <- rep(NA, nalpha)
    for (i in seq(nalpha)) {
      out[i] <- sum(x <= alphaGrid[i])/ntests
    }
    out
  }
  lapply(x, EmpiricalPower, alphaGrid)
}

PowerPlot <- function(ep, alphaGrid, main = "", cex.lab = 1, cex.main = 1, 
                      cex.axis = 1, martext = "") {
  plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1), main = main,
       ylab = "power", xlab = expression(alpha), cex.lab = cex.lab,
       cex.main = cex.main, cex.axis = cex.axis)
  lines(alphaGrid, ep[[1]], col = 1)
  lines(alphaGrid, ep[[2]], col = 4)
  lines(alphaGrid, ep[[3]], col = 5)
  lines(alphaGrid, ep[[4]], col = 6)
  legend("bottomright", legend = c("randomForest", "extraTrees", "UI-Wilcoxon", "UI-t-test"), 
         text.col = c(1, 4, 5, 6), bty = "n")
  mtext(martext, side = 3, at = 0.01, line = 0.5)
}


#######################################
## Auxiliary functions for Figure 8
#######################################

BoxplotAcrossParticipants <- function(aucs, jit, myboxwex = 0.6, mycex = 0.1,
                                      legendPos = "topleft", ylab = "", 
                                      ylim = c(0, 1), cex.lab = 1,
                                      main = "", mymar = c(5, 4, 4, 2) + 0.1,
                                      oParticipants, cex.legend = 1) {
  par(mar = mymar)
  nAlgos <- dim(aucs)[1]
  nParticipants <- nrow(aucs[1,,])
  dat <- t(aucs[1,,])
  dat <- dat[, match(oParticipants, colnames(dat))]
  boxplot(dat, las = 2, border = 1, at = c(1:nParticipants), 
          pars = list(boxwex = myboxwex), cex = mycex, ylab = ylab,
          main = main, ylim = ylim, xaxt = "n", cex.lab = cex.lab,
          outline = FALSE)
  for (i in seq(2, nAlgos, by = 1)) {
    dat <- t(aucs[i,,])
    dat <- dat[, match(oParticipants, colnames(dat))]
    boxplot(dat, add = TRUE, border = i, las = 2, 
            at = c(1:nParticipants) + (1 + i/nAlgos) * jit * (-1)^i, 
            pars = list(boxwex = myboxwex), 
            xaxt = "n", cex = mycex, outline = FALSE)  
  }
  dat <- t(aucs[1,,])
  dat <- dat[, match(oParticipants, colnames(dat))]
  boxplot(dat, add = TRUE, las = 2, border = 1, at = c(1:nParticipants), 
          pars = list(boxwex = myboxwex), cex = mycex,
          main = main, xaxt = "n", outline = FALSE)
  legend(legendPos, legend = dimnames(aucs)[[1]], bty = "n", 
         text.col = seq(nAlgos), cex = cex.legend)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

GetPvaluesAcrossAlgosAndParticipants <- function(x, mtmethod = "BH") {
  nAlgo <- length(x)
  nParticipant <- nrow(x[[1]][[1]])
  pvals <- matrix(NA, nParticipant, nAlgo)
  dimnames(pvals) <- list(names(x[[1]][[2]]), names(x))
  for (i in seq(nParticipant)) {
    for (j in seq(nAlgo)) {
      pvals[i, j] <- x[[j]][[2]][[i]]$pval
    }
  }
  apvals <- apply(pvals, 2, p.adjust, mtmethod)
  list(pvals = pvals, apvals = apvals)
}

PvaluePlotAcrossParticipants <- function(x, 
                                         jit, 
                                         legendPos = "topright", 
                                         ylab = "", 
                                         main = "", 
                                         mymar = c(5, 4, 4, 2) + 0.1, 
                                         oParticipants,
                                         cex.participant = 1,
                                         line = 0.75,
                                         cex = 1,
                                         cex.lab = 1,
                                         cex.legend = 1) {
  par(mar = mymar)
  dat <- -log(x, 10)
  dat <- dat[match(oParticipants, rownames(dat)),]
  ylim <- c(min(dat), max(dat))
  plot(dat[, 1], type = "n", 
       ylab = expression(paste(-log[10], "(p-value)", sep = " ")), 
       xlab = "", xaxt = "n", ylim = ylim, cex.lab = cex.lab)
  pos <- seq(nrow(dat))
  pos2 <- seq(ncol(dat))
  axis(side = 1, at = pos, labels = FALSE)
  mtext(rownames(dat), side = 1, at = seq(nrow(dat)), las = 2, line = line,
        cex = cex.participant)
  for (i in pos) {
    abline(v = pos[i], col = "grey", lty = 3)
    for (j in pos2) {
      points(i + jit * (-1)^j, dat[i, j], col = j, pch = 20, cex = cex)
    }
  }
  legend(legendPos, legend = colnames(x), bty = "n", text.col = seq(pos2),
         cex = cex.legend)
}

GetAurocAcrossAlgosAndParticipants <- function(x) {
  nAlgo <- length(x)
  nParticipant <- nrow(x[[1]][[1]])
  auroc <- matrix(NA, nParticipant, nAlgo)
  dimnames(auroc) <- list(names(x[[1]][[2]]), names(x))
  for (i in seq(nParticipant)) {
    for (j in seq(nAlgo)) {
      auroc[i, j] <- x[[j]][[2]][[i]]$auc
    }
  }
  auroc
}

AurocPlotAcrossParticipants <- function(x, 
                                        jit, 
                                        legendPos = "topright", 
                                        ylab = "", 
                                        main = "", 
                                        mymar = c(5, 4, 4, 2) + 0.1, 
                                        oParticipants,
                                        cex.participant = 1,
                                        line = 0.75,
                                        cex = 1,
                                        cex.lab = 1,
                                        cex.legend = 1) {
  par(mar = mymar)
  dat <- x
  dat <- dat[match(oParticipants, rownames(dat)),]
  ylim <- c(min(dat), max(dat))
  plot(dat[, 1], type = "n", 
       ylab = ylab, 
       xlab = "", xaxt = "n", ylim = ylim, cex.lab = cex.lab)
  pos <- seq(nrow(dat))
  axis(side = 1, at = pos, labels = FALSE)
  for (i in pos) {
    abline(v = pos[i], col = "grey", lty = 3)
    for (j in 2:ncol(dat)) {
      points(i + jit * (-1)^j, dat[i, j], col = j, pch = 20, cex = cex)
    }
    points(i + jit * (-1)^1, dat[i, 1], col = 1, pch = 20, cex = cex)
  }
  legend(legendPos, legend = colnames(x), bty = "n", text.col = seq(ncol(dat)),
         cex = cex.legend)
}

CountBeforeAfterMedicationPerParticipant <- function(x, beforeThr, afterThr) {
  participantIds <- unique(x$healthCode)
  nParticipants <- length(participantIds)
  counts <- data.frame(matrix(NA, nParticipants, 7))
  colnames(counts) <- c("healthCode", "nBefore", "nAfter", "nOther", "nDon't", "nNAs", "n")
  counts[, 1] <- participantIds
  for (i in seq(nParticipants)) {
    pdat <- x[which(x$healthCode == participantIds[i]),]
    aux <- as.character(pdat$momentInDayFormat.json.choiceAnswers)
    counts[i, "n"] <- nrow(pdat)
    counts[i, "nBefore"] <- sum(aux == "Immediately before Parkinson medication")
    counts[i, "nAfter"] <- sum(aux == "Just after Parkinson medication (at your best)")
    counts[i, "nOther"] <- sum(aux == "Another time")
    counts[i, "nDon't"] <- sum(aux == "I don't take Parkinson medications")
    counts[i, "nNAs"] <- counts[i, 7] - sum(counts[i, 2:5])
  }
  idx <- which(counts$nBefore >= beforeThr & counts$nAfter >= afterThr)
  sel <- counts$healthCode[idx]
  nsel <- length(sel)
  nba <- matrix(NA, nsel, 2)
  rownames(nba) <- sel
  colnames(nba) <- c("n_b", "n_a")
  for (i in seq(nsel)) {
    nba[i, 1] <- counts[which(counts$healthCode == sel[i]), "nBefore"]
    nba[i, 2] <- counts[which(counts$healthCode == sel[i]), "nAfter"]
  }
  list(counts = counts, sel = sel, nba = nba) 
}


#######################################
## Auxiliary functions for Figure 9
#######################################

GetParticipantImp <- function(x, id, featnms, top = 3) {
  idx <- match(id, names(x))
  xx <- x[[idx]]
  imp <- xx$imp
  imp <- sort(imp, decreasing = TRUE)
  topfeats <- names(imp)[1:top]
  list(imp = imp, topfeats = topfeats)
}

GetUiTestTopFeatures <- function(x, id, type = c("ttest", "wilcox"), top = 3) {
  idx <- match(id, names(x))
  xx <- x[[idx]]
  pvals <- switch(type, ttest = xx$pvalsTtest, wilcox = xx$pvalsWilcox)
  pvals <- sort(pvals)
  topfeats <- names(pvals)[1:top]
  list(pvals = pvals, topfeats = topfeats)
}

