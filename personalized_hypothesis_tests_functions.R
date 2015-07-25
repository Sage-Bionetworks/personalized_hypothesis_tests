###########################################################
## Functions and utility functions implementing the 
## personalized hypothesis tests.
##
## Ref:
## Chaibub Neto et al. (2015) Personalized hypothesis 
## tests for detecting medication response in Parkinson's
## patients using iPhone sensor data. 
###########################################################

##########################################
## Load required R packages.
##########################################

library(ROCR)
library(randomForest)
library(glmnet)
library(extraTrees)
library(synapseClient)


##########################################
## Utility functions for shaping the
## data in the right format for the  
## personalized tests functions.
##########################################

GetParticipantBeforeAfterData <- function(x, 
                                          participantId) {
  #################################################
  ## Gets the data collected before and after
  ## medication, ignoring data collected at 
  ## another time, or with NAs for the medication
  ## labels.
  ##
  ## Inputs:
  ## x: data.frame containing medication labels, 
  ##    extracted features and other variables
  ## participantIds: character vector with 
  ##                 participant ids (healthCode)
  ##
  ## Output:
  ## pdat: data.frame with a participant 
  ##       before/after medication data
  ##################################################
  pdat <- x[which(x$healthCode == participantId),]
  aux <- as.character(pdat$momentInDayFormat.json.choiceAnswers)
  iBefore <- which(aux == "Immediately before Parkinson medication")
  iAfter <- which(aux == "Just after Parkinson medication (at your best)")
  ii <- sort(c(iBefore, iAfter))
  pdat <- pdat[ii,]
  aux <- as.character(pdat$momentInDayFormat.json.choiceAnswers)
  iBefore <- which(aux == "Immediately before Parkinson medication")
  iAfter <- which(aux == "Just after Parkinson medication (at your best)")
  aux[iBefore] <- "before"
  aux[iAfter] <- "after"  
  pdat$momentInDayFormat.json.choiceAnswers <- factor(aux)
  pdat <- pdat[order(pdat$createdOn),]
  pdat
}

TrainTestRandomSplit <- function(dat, 
                                 nSplits = 2) {
  ##################################################
  ## Split the data into training and testing sets.
  ## 
  ## Inputs:
  ## dat: data matrix or data.frame
  ## nSplits: number of splits (if nSplits = k, 
  ##          approximately 1 fold of the data
  ##          is assigned to the test set, and
  ##          the remaining k-1 folds are used
  ##          for training)
  ##
  ## Outputs:
  ## trainDat: training data
  ## testDat: testing data
  ###################################################
  n <- nrow(dat)
  itest <- sort(sample(n, round(n/nSplits), replace = FALSE))
  itrain <- seq(n)[-itest]
  list(trainDat = dat[itrain,], testDat = dat[itest,])
}

RemoveNAs <- function(x) {
  #####################################
  ## Remove data rows missing one or
  ## more entries.
  ##
  ## Input:
  ## x: matrix or data.frame with 
  ##    medication labels and features
  ##
  ## Output:
  ## x: matrix or data.frame with
  ##    no missing data
  ######################################
  aux <- apply(is.na(x), 1, sum)
  if (sum(aux) != 0) {
    x <- x[-which(aux > 0),]
  }
  x
}


##############################################
## Functions for computing the p-values for
## H_0: AUROC = 0.5 versus H_1: AUROC > 0.5.
##############################################

PvalAUC <- function(auc, 
                    ytest, 
                    predProbs) {
  ################################################
  ## Computes the exact p-value when there are
  ## no ties in the predicted class probabilities,
  ## or the approximated p-value in the presence
  ## of ties.
  ##
  ## Inputs:
  ## auc: observed AUROC
  ## ytest: true labels from the test set
  ##        (used to determine n_b and n_a)
  ## predProbs: predicted class probabilities
  ##
  ## Output:
  ## exact or approximated p-value
  ################################################
  ytest <- factor(ytest)
  ylevels <- levels(ytest)
  pval <- NULL
  if (length(ylevels) == 2) {
    n1 <- sum(ytest == ylevels[1])
    n2 <- sum(ytest == ylevels[2])
  } 
  u <- unique(predProbs)
  if (length(u) == length(predProbs)) { ## get exact p-value
    U <- n1 * n2 * (1 - auc)
    pval <- pwilcox(U, n1, n2, lower.tail = TRUE)
  }
  else { ## get approximated p-value
    ties <- GetTieStats(predProbs)
    n <- n1 + n2
    m <- 1/2
    v <- (n + 1)/(12 * n1 * n2) - ties[[2]]/(12 * n1 * n2 * n * (n - 1))
    pval <- pnorm(auc, m, sqrt(v), lower.tail = FALSE)
  }
  pval
}

GetTieStats <- function(x) {
  ################################################
  ## Computes the variance correction term used 
  ## in the asymptotic normal approximation 
  ## (eq. 4 on manuscript) for the Wilcoxon
  ## rank sum test in the presence of ties.
  ##
  ## Input:
  ## x: predicted class probabilities
  ##
  ## Outputs:
  ## tj: vector with number of ties on group j
  ## aux: correction term
  #################################################
  u <- unique(x)
  idupli <- which(duplicated(x))
  ud <- unique(x[idupli])
  tau <- length(ud)
  tj <- rep(NA, tau)
  for (i in seq(tau)) {
    tj[i] <- sum(x == ud[i])
  }
  list(tj = tj, aux = sum(tj * (tj - 1) * (tj + 1)))
}


###########################################################
## Functions implementing the personalized
## classifier tests.
###########################################################
##
## The following input and output arguments are common 
## to all functions implementing  personalized classifier 
## tests:
##
## Inputs:
##
## nRuns: number of times that we run the classifier 
##        (since we need to split the data into training
##        and test sets, we do this several times and 
##        use the median of the results in order to obtain
##        more robust results with respect to train/test 
##        splits
##
## dat: data.frame with the data from a particular
##      participant
##
## nSplits: number of splits (if nSplits = k, 
##          approximately 1 fold of the data
##          is assigned to the test set, and
##          the remaining k-1 folds are used
##          for training)
##
## splitSeeds: vector of random seeds (with length nRuns)
##             used to that all classifiers are using
##             the same train/test data split at each run
##
## respName: character string with the name of the outcome 
##           variable
##
## refClass: name of the positive class (based on R's 
##           built-in < relation, where the negative and 
##           positive classes are given, respectively, by
##           the smaller and larger "values". For instance,
##           0 < 1, "a" < "b", "after" < "before", and
##           FALSE < TRUE
##
## Outputs:
##
## aucs: vector with AUROCs (one from each of the nRuns)
##
## pvals: vector with p-values (one from each of the nRuns)
##
## auc: median of aucs
##
## pval: median of pvals
############################################################

RandomForestTest <- function(nRuns, 
                             dat, 
                             nSplits = 2, 
                             splitSeeds = NULL, 
                             respName, 
                             refClass = "before") {
  #########################################################
  ## Implements the personalized classifier test based on
  ## the random forest classifier.
  ##
  ## Inputs: see above
  ##
  ## Outputs:
  ## aucs: see above
  ## pvals: see above
  ## Imp: p by nRuns matrix with the importance scores,
  ##      where p represents the number features, and 
  ##      each column stores the importance values from
  ##      one run
  ## auc: see above
  ## pval: see above
  ## imp: mean importance across all nRuns
  ##########################################################
  respColumn <- which(colnames(dat) == respName)
  aucs <- rep(NA, nRuns)
  pvals <- rep(NA, nRuns)
  Imp <- matrix(NA, ncol(dat) - 1, nRuns)
  rownames(Imp) <- colnames(dat)[-respColumn]
  for (i in seq(nRuns)) {
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplit(dat, nSplits)
    trainDat <- RemoveNAs(splitDat$trainDat)
    testDat <- RemoveNAs(splitDat$testDat)
    myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
    fit <- randomForest(myFormula, data = trainDat)
    Imp[, i] <- fit$importance[, 1]
    pred <- predict(fit, testDat[, -respColumn], type = "prob")[, refClass]
    ytest <- testDat[, respColumn]
    predobj <- prediction(pred, ytest)
    aucs[i] <- performance(predobj, "auc")@y.values[[1]]
    pvals[i] <- PvalAUC(aucs[i], ytest, as.numeric(pred))    
  }
  mauc <- median(aucs)
  pval <- median(pvals)
  imp <- apply(Imp, 1, mean)
  list(pval = pval, auc = mauc, imp = imp, aucs = aucs, pvals = pvals, Imp = Imp)
}

ExtraTreesTest <- function(nRuns, 
                           dat, 
                           nSplits = 2, 
                           splitSeeds, 
                           respName, 
                           refClass = "before") {
  #########################################################
  ## Implements the personalized classifier test based on
  ## the extra tree classifier.
  ##
  ## Inputs: see above
  ##
  ## Outputs: see above
  ##########################################################
  respColumn <- which(colnames(dat) == respName)
  aucs <- rep(NA, nRuns)
  pvals <- rep(NA, nRuns)
  for (i in seq(nRuns)) {
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplit(dat, nSplits)
    trainDat <- RemoveNAs(splitDat$trainDat)
    testDat <- RemoveNAs(splitDat$testDat)
    myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
    fit <- extraTrees(x = trainDat[, -respColumn], y = trainDat[, respColumn])   
    pred <- predict(fit, testDat[, -respColumn], probability = TRUE)[, refClass]
    ytest <- testDat[, respColumn]
    predobj <- prediction(pred, ytest)
    aucs[i] <- performance(predobj, "auc")@y.values[[1]]
    pvals[i] <- PvalAUC(aucs[i], ytest, as.numeric(pred))    
  }
  mauc <- median(aucs)
  pval <- median(pvals)
  list(pval = pval, auc = mauc, aucs = aucs, pvals = pvals)
}

GlmTest <- function(nRuns, 
                    dat,  
                    nSplits = 2, 
                    splitSeeds,
                    respName) {
  #########################################################
  ## Implements the personalized classifier test based on
  ## the logist regression classifier.
  ##
  ## Inputs: see above
  ##
  ## Outputs:
  ## aucs: see above
  ## pvals: see above
  ## betas: p by nRuns matrix with regr. coef. estimates,
  ##        where p represents the number features, and 
  ##        each column stores the coefficients from
  ##        one run
  ## auc: see above
  ## pval: see above
  ## beta: mean regr. coef. estimates across all nRuns
  ##########################################################
  aucs <- rep(NA, nRuns)
  pvals <- rep(NA, nRuns)
  respColumn <- which(colnames(dat) == respName)
  betas <- matrix(NA, ncol(dat) - 1, nRuns)
  rownames(betas) <- colnames(dat)[-respColumn]
  myFormula <- as.formula(paste(respName, " ~ .", sep = ""))
  for (i in seq(nRuns)) {
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplit(dat, nSplits)
    trainDat <- RemoveNAs(splitDat$trainDat)
    testDat <- RemoveNAs(splitDat$testDat)
    ytest <- testDat[, respColumn] 
    Xtest <- testDat[, -respColumn]
    fit <- glm(myFormula, data = trainDat, family = "binomial")
    pred <- predict(fit, newdata = testDat[, -respColumn], type = "response")
    predobj <- prediction(pred, ytest)
    aucs[i] <- performance(predobj, "auc")@y.values[[1]]
    pvals[i] <- PvalAUC(aucs[i], ytest, as.numeric(pred)) 
    betas[, i] <- fit$coef[-1]
  }  
  mauc <- median(aucs)
  pval <- median(pvals)
  mbeta <- apply(betas, 1, mean)
  list(pval = pval, auc = mauc, beta = mbeta, aucs = aucs, pvals = pvals, betas = betas)
}

EnetTest <- function(nRuns, 
                     dat, 
                     nSplits = 2, 
                     splitSeeds,
                     respName,
                     alphaGrid = seq(0.1, 0.9, by = 0.1),
                     nfolds = 3) {
  #########################################################
  ## Implements the personalized classifier test based on
  ## the elastic-net penalized logist regression classifier.
  ##
  ## Inputs:
  ## nRuns, dat, respName, nSplits, splitSeeds: see above
  ## alphaGrid: grid for elatic-net's alpha tuning parameter
  ## nfolds: number of cross-validation folds used for
  ##         tunning parameter optimization
  ##
  ## Outputs:
  ## aucs: see above
  ## pvals: see above
  ## betas: p by nRuns matrix with regr. coef. estimates,
  ##        where p represents the number features, and 
  ##        each column stores the coefficients from
  ##        one run
  ## auc: see above
  ## pval: see above
  ## beta: mean regr. coef. estimates across all nRuns
  ##########################################################  
  EnetClassFit <- function(trainDat, 
                           Xtest,
                           respName,
                           nfolds = 10,
                           alphaGrid = seq(0.1, 0.9, by = ),
                           cvSeed = NULL) {
    if (!is.null(cvSeed)) {
      set.seed(cvSeed)
    }
    respColumn <- which(colnames(trainDat) == respName)
    ytrain <- trainDat[, respColumn]
    Xtrain <- trainDat[, -respColumn]
    featNames <- colnames(Xtrain)
    myform <- as.formula(paste(" ~ -1 + ", paste(featNames, collapse = " + "), sep = ""))
    Xtrain <- model.matrix(myform, Xtrain)
    Xtest <- model.matrix(myform, Xtest)
    nAlpha <- length(alphaGrid)
    lambs <- rep(NA, nAlpha)
    alphaError <- rep(NA, nAlpha)
    for (j in seq(nAlpha)) {
      cvFit <- cv.glmnet(Xtrain, ytrain, nfolds = nfolds, alpha = alphaGrid[j], 
                         family = "binomial", type.measure = "auc")
      lambs[j] <- cvFit$lambda.min
      alphaError[j] <- cvFit$cvm[which(cvFit$lambda == cvFit$lambda.min)]
    }
    best <- which.min(alphaError)
    bestAlpha <- alphaGrid[best]
    bestLambda <- lambs[best]
    bestFit <- glmnet(Xtrain, ytrain, alpha = bestAlpha, lambda = bestLambda, 
                      family = "binomial") 
    pred <- predict(bestFit, newx = Xtest, s = bestLambda, type = "response")
    betas <- bestFit$beta[, 1]
    list(pred = pred, betas = betas)    
  }
  aucs <- rep(NA, nRuns)
  pvals <- rep(NA, nRuns)
  respColumn <- which(colnames(dat) == respName)
  betas <- matrix(NA, ncol(dat) - 1, nRuns)
  rownames(betas) <- colnames(dat)[-respColumn]
  for (i in seq(nRuns)) {
    set.seed(splitSeeds[i])
    splitDat <- TrainTestRandomSplit(dat, nSplits)
    trainDat <- RemoveNAs(splitDat$trainDat)
    testDat <- RemoveNAs(splitDat$testDat)
    ytest <- testDat[, respColumn] 
    Xtest <- testDat[, -respColumn]
    fit <- EnetClassFit(trainDat, Xtest, respName, nfolds = nfolds, 
                        alphaGrid = alphaGrid, cvSeed = NULL)
    predobj <- prediction(fit$pred, ytest)
    aucs[i] <- performance(predobj, "auc")@y.values[[1]]
    pvals[i] <- PvalAUC(aucs[i], ytest, as.numeric(fit$pred)) 
    betas[, i] <- fit$betas
  }  
  mauc <- median(aucs)
  pval <- median(pvals)
  mbeta <- apply(betas, 1, mean)
  list(pval = pval, auc = mauc, beta = mbeta, aucs = aucs, pvals = pvals, betas = betas)
}

RunPerClassTests <- function(selParticipants,
                             nRuns, 
                             dat, 
                             nSplits,
                             splitSeeds,
                             respName, 
                             featNames,
                             ClassTest, ...) {
  ##########################################################
  ## Runs the personalized classification tests across a
  ## set of selected participants.
  ##
  ## Inputs:
  ## selParticipants: character vector with the participant 
  ##                  ids (healthCode)
  ## dat: data.frame containing the data from all participants
  ## nRuns, nSplits, splitSeeds, respName: see above
  ## featNames: vector with the names of the features to be
  ##            included in the analysis
  ## ClassTest: name of the function implementing the
  ##            personalized classifier test (e.g., 
  ##            RandomForestTest, ExtraTreesTest,
  ##            GlmTest, EnetTest)
  ## ... : additional parameters required by the ClassTest
  ##       function
  ##
  ## Outputs:
  ## AUCpvals: matrix reporting the AUROC and associated 
  ##           p-value, across all selected participants
  ## participantOutputs: list containing the output of the
  ##                     ClassTest function of each one of
  ##                     the selected participants
  ########################################################### 
  call <- match.call()
  nParticipants <- length(selParticipants)
  participantOutputs <- vector(mode = "list", length = nParticipants)
  names(participantOutputs) <- selParticipants
  AUCpvals <- matrix(NA, nParticipants, 2)
  colnames(AUCpvals) <- c("AUC", "pvalue")
  rownames(AUCpvals) <- selParticipants
  respCol <- match(respName, colnames(dat))
  featCols <- match(featNames, colnames(dat))
  for (i in seq(nParticipants)) {
    cat("running participant ", i, "\n")  
    pdat <- GetParticipantBeforeAfterData(dat, selParticipants[i])[, c(respCol, featCols)]
    aux <- ClassTest(nRuns, pdat, nSplits, splitSeeds, respName, ...)
    participantOutputs[[i]] <- aux
    AUCpvals[i, 1] <- aux$auc
    AUCpvals[i, 2] <- aux$pval
  }
  oo <- order(AUCpvals[, "pvalue"], decreasing = FALSE)
  AUCpvals <- AUCpvals[oo,] 
  list(AUCpvals = AUCpvals, participantOutputs = participantOutputs)
}


##################################
## Functions implementing the 
## Union-Intersection tests
##################################

UITests <- function(dat, 
                    mtmethod = "BH", 
                    alternative = "two.sided", 
                    respName) {
  ###########################################################
  ## Implements the UI-tests based on t-tests and Wilcoxon
  ## rank sum tests.
  ##
  ## Inputs:
  ## dat: data.frame with the data from a particular
  ##      participant
  ## mtmethod: multiple testing correction method (need 
  ##           to be one of the methods implemented in
  ##           the p.adjust() function - our default,
  ##           BH, corresponds to the Benjamini-Hochberg
  ##           method)
  ## alternative: character string specifying the 
  ##              alternative hypothesis
  ## respName: character string with the name of the 
  ##           outcome variable
  ##
  ## Outputs:
  ## pvalTtest: p-value of the UI-test based on t-tests
  ## pvalWilcox: p-value of the UI-test based on 
  ##             Wilcoxon's rank sum tests
  ## pvalsTtest: t-test p-values for each separate 
  ##             feature
  ## pvalsWilcox: Wilcoxon's test p-values for each
  ##              separate feature
  ###########################################################
  respColumn <- which(colnames(dat) == respName)
  nfeat <- ncol(dat) - 1
  featnms <- colnames(dat)[-respColumn]
  featColumns <- seq(ncol(dat))[-respColumn]
  pvalsTtest <- rep(NA, nfeat)
  pvalsWelchs <- rep(NA, nfeat)
  pvalsWilcox <- rep(NA, nfeat)
  names(pvalsTtest) <- featnms
  names(pvalsWilcox) <- featnms
  respLevels <- levels(dat[, respColumn])
  for (i in featColumns) {
    xa <- dat[, i][which(dat[, respColumn] == respLevels[1])]
    xb <- dat[, i][which(dat[, respColumn] == respLevels[2])] 
    xa <- xa[which(!is.na(xa))]
    xb <- xb[which(!is.na(xb))]
    pvalsTtest[i-1] <- t.test(xa, xb, var.equal = TRUE, alternative = alternative)$p.value
    pvalsWilcox[i-1] <- wilcox.test(xa, xb, alternative = alternative)$p.value
  }
  pvalTtest <- min(p.adjust(pvalsTtest, method = mtmethod))
  pvalWilcox <- min(p.adjust(pvalsWilcox, method = mtmethod))
  list(pvalTtest = pvalTtest, 
       pvalWilcox = pvalWilcox,
       pvalsTtest = pvalsTtest,  
       pvalsWilcox = pvalsWilcox)
}

RunUITests <- function(selParticipants, 
                       dat, 
                       sort.by = "none") {
  ###########################################################
  ## Runs the Union-Intersection tests across a set of 
  ## selected participants.
  ##
  ## Inputs:
  ## selParticipants: character vector with the participant 
  ##                  ids (healthCode)
  ## dat: data.frame containing the data from all 
  ##      participants
  ## sort.by: sort the results by the UI-t-test ("t.test"),
  ##          or by the UI-Wilcoxon test ("Wilcoxon"),
  ##          or leave results unsorted ("none")
  ##
  ## Outputs:
  ## UIpvals: matrix reporting the p-values for the UI-t-test 
  ##          and UI-Wilcoxon test, across all selected 
  ##          participants
  ## participantOutputs: list containing the output of the
  ##                     UItests function of each one of
  ##                     the selected participants
  ###########################################################
  nParticipants <- length(selParticipants)
  participantOutputs <- vector(mode = "list", length = nParticipants)
  names(participantOutputs) <- selParticipants
  UIpvals <- matrix(NA, nParticipants, 2)
  colnames(UIpvals) <- c("t.test", "Wilcoxon")
  rownames(UIpvals) <- selParticipants
  respName <- "momentInDayFormat.json.choiceAnswers"
  for (i in seq(nParticipants)) {
    cat("running participant ", i, "\n") 
    pdat <- GetParticipantBeforeAfterData(dat, selParticipants[i])[, 8:32]
    aux <- UITests(pdat, mtmethod = "BH", alternative = "two.sided", respName)
    participantOutputs[[i]] <- aux
    UIpvals[i, 1] <- aux$pvalTtest
    UIpvals[i, 2] <- aux$pvalWilcox
  }
  if (sort.by != "none") {
    oo <- order(UIpvals[, sort.by], decreasing = FALSE)
    UIpvals <- UIpvals[oo,] 
  }
  list(UIpvals = UIpvals, participantOutputs = participantOutputs)
}


####################################################
## Utility functions used for simulating data and
## running the statistical power simulation study.
####################################################

SimMultipleSeasonalSeries <- function(n, 
                                      a, 
                                      sigs, 
                                      err = 0.1) {
  t <- seq(n)
  nSeries <- length(sigs)
  x <- a * cos(pi*t)
  y <- rep(NA, n)
  y[which(x > 0)] <- "a"
  y[which(x < 0)] <- "b"
  idx1 <- which(x > 0)
  y[sample(idx1, round(n * err/2), replace = FALSE)] <- "b"
  idx2 <- which(x < 0)
  y[sample(idx2, round(n * err/2), replace = FALSE)] <- "a"    
  y <- factor(y)
  X <- matrix(NA, n, nSeries)
  colnames(X) <- paste("x", seq(nSeries), sep = "")
  for (j in seq(nSeries)) {
    X[, j] <- x + rnorm(n, 0, sigs[j])
  }
  data.frame(y, X)
}

RunSimulations <- function(nSim, 
                           n, 
                           p, 
                           a, 
                           maxsig, 
                           pmaxsig, 
                           err, 
                           myseed, 
                           nSplits = 2) {
  pvalRF <- rep(NA, nSim)
  pvalET <- rep(NA, nSim)
  pvalW <- rep(NA, nSim)
  pvalT <- rep(NA, nSim)
  set.seed(myseed)
  datSeeds <- sample(10000:100000, nSim, replace = FALSE)
  for (i in seq(nSim)) {
    cat(i, "\n")
    set.seed(datSeeds[i])
    sigs <- c(seq(1, maxsig, length.out = pmaxsig), rep (maxsig, p - pmaxsig))
    dat <- SimMultipleSeasonalSeries(n, a, sigs, err = err)
    ui <- UITests(dat, mtmethod = "BH", alternative = "two.sided", respName = "y")
    rf <- RandomForestTest(nRuns = 1, dat, nSplits, splitSeed = NULL, respName = "y", refClass = "b")
    et <- ExtraTreesTest(nRuns = 1, dat, nSplits, splitSeed = NULL, respName = "y", refClass = "b")
    pvalRF[i] <- rf$pval
    pvalET[i] <- et$pval
    pvalW[i] <- ui$pvalWilcox
    pvalT[i] <- ui$pvalTtest
  }
  list(pvalRF = pvalRF, pvalET = pvalET, pvalW = pvalW, pvalT = pvalT)
}
