############################################################
## Code used for generating Figures 1 to 9 on the paper:
##
## Chaibub Neto et al. (2015) Personalized hypothesis 
## tests for detection medication response in Parkinson's
## patients using iPhone sensor data. 
############################################################

library(devtools)

## source auxiliary functions for shaping and plotting the 
## analysis results
##
#source("personalized_hypothesis_tests_auxiliary_functions.R")
source_url("https://raw.githubusercontent.com/Sage-Bionetworks/personalized_hypothesis_tests/master/personalized_hypothesis_tests_auxiliary_functions.R?token=ABoVeWCp1hIA9WssBCJEbtD-M2G5lx_7ks5VvLaDwA%3D%3D")

## source functions for shaping the data and running
## the personalized hypotheis tests
##
#source("personalized_hypothesis_tests_functions.R")
source_url("https://raw.githubusercontent.com/Sage-Bionetworks/personalized_hypothesis_tests/master/personalized_hypothesis_tests_functions.R?token=ABoVeay25neuXIVrxREoWsVq7bxjautrks5VvLRawA%3D%3D")

## load the tapping data on the subset of 57 patients
## that performed at least 30 tapping tasks before medication
## and 30 tapping tasks after medication
##
#load("patients_subset_tapping_6_25_15.RData")
load(synGet("syn4649829")@filePath) ## or load directly from Synapse

## path to the folder storing the figures
##
fig.path <- ""

########################################
## Figure 1
########################################

library(MASS)

## generate 10,000 simulated null data sets from
## model (i), and apply a t-test to each one of them
##
nSim <- 10000
n1 <- 30
n2 <- 30
mu1 <- mu2 <- 0
rho1 <- 0
rho2 <- 0
Sig1 <- matrix(rho1, n1, n1)
diag(Sig1) <- 1
Sig2 <- matrix(rho2, n2, n2)
diag(Sig2) <- 1
X1 <- matrix(NA, nSim, n1 + n2)
tpvals1 <- rep(NA, nSim)
set.seed(987654321)
for (i in seq(nSim)) {
  cat(i, "\n")
  x1 <- mvrnorm(1, rep(mu1, n1), Sig1)
  x2 <- mvrnorm(1, rep(mu2, n2), Sig2)
  tpvals1[i] <- t.test(x1, x2, var.equal = TRUE)$p.value
  X1[i,] <- c(x1, x2)
}

## generate 10,000 simulated null data sets from
## model (ii), and apply a t-test to each one of them 
##
rho1 <- 0.95
rho2 <- 0.95
Sig1 <- matrix(rho1, n1, n1)
diag(Sig1) <- 1
Sig2 <- matrix(rho2, n2, n2)
diag(Sig2) <- 1
X2 <- matrix(NA, nSim, n1 + n2)
tpvals2 <- rep(NA, nSim)
set.seed(12345)
for (i in seq(nSim)) {
  cat(i, "\n")
  x1 <- mvrnorm(1, rep(mu1, n1), Sig1)
  x2 <- mvrnorm(1, rep(mu2, n2), Sig2)
  tpvals2[i] <- t.test(x1, x2, var.equal = TRUE)$p.value
  X2[i,] <- c(x1, x2)
}

## generate 10,000 simulated null data sets from
## model (iii), and apply a t-test to each one of them 
##
n <- n1 + n2
mu <- 0
rho <- 0.95
Sig <- matrix(rho, n, n)
diag(Sig) <- 1
i1 <- seq(n1)
i2 <- seq(n1 + 1, n, by = 1)
X3 <- matrix(NA, nSim, n)
tpvals3 <- rep(NA, nSim)
set.seed(12345)
for (i in seq(nSim)) {
  cat(i, "\n")
  x <- mvrnorm(1, rep(mu, n), Sig)
  tpvals3[i] <- t.test(x[i1], x[i2], var.equal = TRUE)$p.value
  X3[i,] <- x
}

## generate Figure 1
##
cl <- 1.2; ca <- 1.2; cm <- 1.4
k <- 30 ## show top 30 simulations only
postscript(paste(fig.path, "fig1.ps", sep = ""), width = 10, height = 4)
par(mfrow = c(2, 3), mar = c(1.5, 3.5, 1.5, 0.5) + 0.1, mgp = c(2, 0.75, 0))
image(t(X1[1:k,]), xaxt = "n", yaxt = "n", main = "independ. within and between groups", cex.main = cm)
abline(v = 0.5)
mtext("(a)", side = 3, at = 0.9, cex = 1.2, line = -2)
mtext("group 1", side = 1, at = 0.25, line = 1)
mtext("group 2", side = 1, at = 0.75, line = 1)
mtext("simulated data sets", side = 2, at = 0.5, line = 1)
image(t(X2[1:k,]), xaxt = "n", yaxt = "n", main = "dep. within and indep. between groups", cex.main = cm)
abline(v = 0.5)
mtext("(b)", side = 3, at = 0.9, cex = 1.2, line = -2)
mtext("group 1", side = 1, at = 0.25, line = 1)
mtext("group 2", side = 1, at = 0.75, line = 1)
mtext("simulated data sets", side = 2, at = 0.5, line = 1)
image(t(X3[1:k,]), xaxt = "n", yaxt = "n", main = "depend. within and between groups", cex.main = cm)
abline(v = 0.5)
mtext("(c)", side = 3, at = 0.9, cex = 1.2, line = -2)
mtext("group 1", side = 1, at = 0.25, line = 1)
mtext("group 2", side = 1, at = 0.75, line = 1)
mtext("simulated data sets", side = 2, at = 0.5, line = 1)
par(mar = c(3, 3.5, 1.5, 0.5) + 0.1)
hist(tpvals1, probability = TRUE, ylim = c(0, 1.4), xlab = "p-values", 
     cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(d)", side = 3, at = 0.9, cex = 1.2, line = -2)
hist(tpvals2, probability = TRUE, xlab = "p-values",
     cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(e)", side = 3, at = 0.9, cex = 1.2, line = -2)
hist(tpvals3, probability = TRUE, ylim = c(0, 1.4), xlab = "p-values",
     cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(f)", side = 3, at = 0.9, cex = 1.2, line = -2)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
dev.off()


########################################
## Figure 2
########################################

## simulate data from a AR(1) process
##
n <- 100
set.seed(123456789)
dat <- arima.sim(list(order = c(1,1,0), ar = -0.5), n = n)
x <- as.vector(dat)[-1]
n1 <- round(n/2)
n2 <- n - n1
set.seed(12345) ## randomly assign the medication labels indexes
before <- sort(sample(seq(n), n1, replace = FALSE)) 
after <- setdiff(seq(n), before)

## generate Figure 2
##
xaxis <- seq(n)
cx <- 1; cl <- 1.5; ca <- 1.2; cm <- 1.4; cx <- 0.8
postscript(paste(fig.path, "fig2.ps", sep = ""), width = 10, height = 2)
par(mfrow = c(1, 3), mar = c(3.4, 3.5, 1.5, 0.5) + 0.1, mgp = c(2, 0.75, 0))
plot(x, type = "l", main = "(a) AR(1) data simulated under the null", 
     ylab = "feature", xlab = "time", cex.main = cm, cex.lab = cl, cex.axis = ca)
points(xaxis[before], x[before], col = "red", pch = 20, cex = cx)
points(xaxis[after], x[after], col = "blue", pch = 20, cex = cx)
legend("bottomright", legend = c("before", "after"), text.col = c("red", "blue"), bty = "n")
acf(x, main = "", cex.main = cm, cex.lab = cl, cex.axis = ca)
mtext("(b) autocorrelation plot of feature data", side = 3, at = 10, cex = 0.9, 
      font = 2, line = 0.2)
plot(x, col = "grey", lwd = 2, type = "l", pch = 20, ylab = "feature", xlab = "time",
     main = "(c) before and after medication series",
     cex.main = cm, cex.lab = cl, cex.axis = ca)
legend("bottomright", legend = c("before", "after"), text.col = c("red", "blue"), bty = "n")
points(xaxis[before], x[before], col = "red", pch = 20, cex = cx)
points(xaxis[after], x[after], col = "blue", pch = 20, cex = cx)
lines(xaxis[before], x[before], col = "red", lwd = 1)
lines(xaxis[after], x[after], col = "blue", lwd = 1)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
dev.off()


########################################
## Figure 3
########################################

## get data from patient cdb18d7b-6793-48c0-8085-63c2a71ce5a9
## collected during April 2015
##
id <- "cdb18d7b-6793-48c0-8085-63c2a71ce5a9"
pdat <- scases[scases$healthCode == id,]
pdat <- pdat[order(pdat$createdOn),]
pdat <- pdat[20:110,]
## generate one shuffled version of the medication labels
sdat <- ShuffleBeforeAfterLabels(pdat, seed = 1235) 

## get data on the number of taps feature before and after medication
##
xb <- pdat[which(pdat$momentInDayFormat.json.choiceAnswers == "Immediately before Parkinson medication"), "numberTaps"]
xa <- pdat[which(pdat$momentInDayFormat.json.choiceAnswers == "Just after Parkinson medication (at your best)"), "numberTaps"]
nb <- length(xb)
na <- length(xa)

## run permutation tests based on the t-test statistic, and 
## get the analytical null distribution
##
set.seed(7654321)
permt <- TTestPerm(B = 1e+4, xa, xb)
xaxist <- seq(min(permt$stat), max(permt$stat), length.out = 100)
tdensi <- dt(xaxist, df = nb + na - 2)

## run permutation tests based on the Wilcoxon rank sum test 
## statistic, and get the (approximate) analytical null 
## distribution (eq. 3 in the paper shows a slightly more 
## accurate approximation, which is used in our classifier
## tests)
##
set.seed(7654321)
permw <- WilcoxTestPerm(B = 1e+4, xa, xb)
xaxisw <- seq(min(permw$stat), max(permw$stat), length.out = 100)
wdensi <- dnorm(xaxisw, na*nb/2, sqrt(na*nb*(na+nb+1)/12))

## get the analytical p-values (t- and Wilcoxon tests)
## from 10,000 null data sets
##
set.seed(1234567)
null <- GetNullModelAnalyticalPvalues(nSim = 1e+4, pdat)

## generate Figure 3
##
cl <- 1.2; ca <- 1.2; cm <- 1.2
postscript(paste(fig.path, "fig3.ps", sep = ""), width = 10, height = 6)
par(mfrow = c(3, 2), mar = c(3.4, 3.5, 1.5, 0.5) + 0.1, mgp = c(1.75, 0.5, 0))
PlotBeforeAfterMedication(healthCodeId = id, featName = "numberTaps", 
                          dat = pdat, ylab = "number of taps", ylim = NULL, 
                          cex.dot = 1.5, legendPos = "topleft",
                          main = "(a)   original data", cex.lab = cl,
                          cex.axis = ca, cex.main = cm)
PlotBeforeAfterMedication(healthCodeId = id, featName = "numberTaps", 
                          dat = sdat, ylab = "number of taps", ylim = NULL, 
                          cex.dot = 1.5, legendPos = "topleft",
                          main = "(b)    shuffled before/after medication labels",
                          cex.lab = cl, cex.axis = ca, cex.main = cm)
hist(permt$stat, nclass = 40, probability = TRUE, 
     main = "permutation null (t-test statistic)", xlab = "test statistic",
     cex.lab = cl, cex.axis = ca)
lines(xaxist, tdensi, col = "red")
mtext("(c)", side = 3, at = -4.5, cex = 1.2, line = -2)
hist(permw$stat, nclass = 30, probability = TRUE, 
     main = "permutation null (Wilcoxon test statistic)", xlab = "test statistic",
     cex.lab = cl, cex.axis = ca)
lines(xaxisw, wdensi, col = "red")
mtext("(d)", side = 3, at = 225, cex = 1.2, line = -2)
hist(null$pvalsT, probability = TRUE, nclass = 20, xlab = "p-values", 
     ylim = c(0, 1.4), cex.lab = cl, cex.axis = ca, cex.main = cm,
     main = "p-value distribution under the null (from analytical t-tests)")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(e)", side = 3, at = 0.05, cex = 1.2, line = -1.5)
hist(null$pvalsW, probability = TRUE, nclass = 20, xlab = "p-values", 
     ylim = c(0, 1.4), cex.lab = cl, cex.axis = ca, cex.main = cm,
     main = "p-value distribution under the null (from analytical Wilcoxon tests)")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(f)", side = 3, at = 0.05, cex = 1.2, line = -1.5)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
dev.off()


########################################
## Figure 4
########################################

id <- "b5d300e5-4397-4f24-ab97-6ae76a3c956b"
pdat <- scases[scases$healthCode == id,]

set.seed(987654321)
o1 <- CheckNullUIPvalDistribution(pdat, id, nperm = 10000, breakFeatCor = FALSE)
set.seed(987654321)
o2 <- CheckNullUIPvalDistribution(pdat, id, nperm = 10000, breakFeatCor = TRUE)

## generate Figure 4
##
nc <- 20
mylim <- c(0, 2.5)
cl <- 1.2; ca <- 1.2; cm <- 1.4
postscript(paste(fig.path, "fig4.ps", sep = ""), width = 10, height = 2)
par(mfrow = c(1, 4), mar = c(3.4, 3.5, 1.5, 0.5) + 0.1, mgp = c(2, 0.75, 0))
hist(o1$pvalsWilcoxM[, 2], probability = TRUE, nclass = nc, ylim = mylim,
     xlab = "p-values", cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "median tapping interval")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(a)", side = 3, at = 0.1, cex = 1.2, line = -2)
hist(o1$pvalsWilcoxM[, 10], probability = TRUE, nclass = nc, ylim = mylim,
     xlab = "p-values", cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "number of taps")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(b)", side = 3, at = 0.1, cex = 1.2, line = -2)
hist(o1$pwilcoxson, probability = TRUE, nclass = nc, ylim = mylim, 
     xlab = "p-values", cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "UI-test with cor. feat.")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(c)", side = 3, at = 0.1, cex = 1.2, line = -2)
hist(o2$pwilcoxson, probability = TRUE, nclass = nc, ylim = mylim, 
     xlab = "p-values", cex.axis = ca, cex.lab = cl, cex.main = cm,
     main = "UI-test with uncor. feat.")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(d)", side = 3, at = 0.1, cex = 1.2, line = -2)
par(mfrow = c(1, 1))
dev.off()


########################################
## Figure 5
########################################

## get before/after medication data for patient
## b5d300e5-4397-4f24-ab97-6ae76a3c956b
##
id <- "b5d300e5-4397-4f24-ab97-6ae76a3c956b"
pdat <- GetParticipantBeforeAfterData(scases, id)[, 8:32]

## generate the permutation null distribution for the 
## AUROC metric (based on the random forest classifier)
##
permRF <- RandomForestPermNull(nperm = 1e+4, pdat, nSplits = 2, splitSeed = 123, 
                               respName = "momentInDayFormat.json.choiceAnswers")
paucs <- permRF$paucs
nb <- permRF$nb
na <- permRF$na
aux <- permRF$aux
n <- nb + na

## get the corresponding (approximated) analytical null
## distribution (eq. 4 in the paper)
##
xaxis <- seq(0, 1, length.out = 1000)
normapprox <- dnorm(xaxis, 1/2, sqrt((n + 1)/(12 * nb * na) - aux/(12 * nb * na * n * (n - 1))))

## get the analytical p-values (from random forestclassifier test) 
## from 10,000 null data sets
##
set.seed(987654321)
null <- GetNullModelAnalyticalPvalues2(nSim = 1e+4, pdat)

## generate Figure 5
##
cl <- 0.9; ca <- 0.9; cm <- 0.8
postscript(paste(fig.path, "fig5.ps", sep = ""), width = 10, height = 2)
par(mfrow = c(1, 2), mar = c(3.4, 3.5, 1.5, 0.5) + 0.1, mgp = c(2, 0.75, 0))
hist(paucs, probability = TRUE, nclass = 20, xlab = "test statistic",
     cex.main = cm, cex.lab = cl, cex.axis = ca,
     main = "permutation null (random forest test)")
lines(xaxis, normapprox, col = "red")
mtext("(a)", side = 3, at = 0.275, cex = 1.2, line = -1)
hist(null, nclass = 20, probability = TRUE, ylim = c(0, 1.4), 
     cex.main = cm, cex.lab = cl, cex.axis = ca, xlab = "p-values",
     main = "p-value distribution under the null (analytic p-values)")
segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, col = "red")
segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, col = "red")
segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, col = "red")
mtext("(b)", side = 3, at = 0.1, cex = 1.2, line = -1)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
dev.off()


########################################
## Figure 6
########################################

N <- 30
t <- seq(N)
A <- 0.25
x <- A * cos(pi*t)
ibefore <- which(x < 0)
iafter <- which(x > 0)
mycex <- 1.5
mylim <- c(-0.8, 0.8)
sig <- 0.3

## generate Figure 6
##
postscript(paste(fig.path, "fig6.ps", sep = ""), width = 10, height = 2)
par(mfrow = c(1, 4), mar = c(3, 2, 2, 0.5) + 0.1, mgp = c(2, 0.75, 0))
plot(x, type = "l", main = "periodic signal", ylim = mylim, 
     xlab = "time", ylab = "")
abline(h = A, col = "grey")
abline(h = -A, col = "grey")
mtext("(a)", side = 3, at = 29, line = -1.5, cex = 1.1)
plot(x, type = "l", main = "add before/after medication labels", 
     ylim = mylim, xlab = "time", ylab = "")
abline(h = A, col = "grey")
abline(h = -A, col = "grey")
points(ibefore, x[ibefore], pch = 20, cex = mycex, col = "red")
points(iafter, x[iafter], pch = 20, cex = mycex, col = "blue")
legend("bottomleft", legend = c("before", "after"), text.col = c("red", "blue"), bty = "n")
mtext("(b)", side = 3, at = 29, line = -1.5, cex = 1.1)
set.seed(123)
plot(x, type = "l", main = "introduce labeling errors", 
     ylim = mylim, xlab = "time", ylab = "")
abline(h = A, col = "grey")
abline(h = -A, col = "grey")
points(ibefore, x[ibefore], pch = 20, cex = mycex, col = "red")
points(iafter, x[iafter], pch = 20, cex = mycex, col = "blue")
legend("bottomleft", legend = c("before", "after"), text.col = c("red", "blue"), bty = "n")
err <- 0.2
i1 <- sample(ibefore, round(N * err/2), replace = FALSE)
points(i1, x[i1], pch = 20, cex = mycex, col = "blue")
i2 <- sample(iafter, round(N * err/2), replace = FALSE)
points(i2, x[i2], pch = 20, cex = mycex, col = "red")
mtext("(c)", side = 3, at = 29, line = -1.5, cex = 1.1)
set.seed(12345)
x2 <- x + rnorm(N, 0, sig) 
plot(x2, type = "l", main = "add random noise", 
     ylim = mylim, xlab = "time", ylab = "")
points(ibefore, x2[ibefore], pch = 20, cex = mycex, col = "red")
points(iafter, x2[iafter], pch = 20, cex = mycex, col = "blue")
points(i1, x2[i1], pch = 20, cex = mycex, col = "blue")
points(i2, x2[i2], pch = 20, cex = mycex, col = "red")
legend("bottomleft", legend = c("before", "after"), text.col = c("red", "blue"), bty = "n")
mtext("(d)", side = 3, at = 29, line = -1.5, cex = 1.1)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
dev.off()


########################################
## Figure 7
########################################

## load results generated by the code 
## on run_power_simulations.R
##
load("power_simulation_study_outputs.RData")
#load(synGet("syn4677520")@filePath) ## or load directly from Synapse

## set the significance level grid
##
alphas <- seq(0, 1, by = 0.001)
nms <- names(powerSim)

## generate Figure 7
##
postscript(paste(fig.path, "fig7.ps", sep = ""), width = 10, height = 4)
par(mfrow = c(2, 3), mar = c(3, 3, 1.5, 1) + 0.1, mgp = c(2, 0.75, 0))
martext <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
for (k in seq(6)) {
  ep <- GetEmpiricalPower(powerSim[[k]], alphas)
  PowerPlot(ep, alphas, main = nms[k], cex.lab = 1.5, 
            cex.main = 1.5, cex.axis = 1.2, martext[k]) 
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(2, 1, 0))
dev.off()


########################################
## Figure 8
########################################

## load results generated by the code 
## on run_real_data_analysis_tapping_6_25_15.R
##
load("real_data_analysis_outputs.RData")
#load(synGet("syn4677553")@filePath) ## or load directly from Synapse

## get participants that performed at least 30 tapping
## tasks before medication, and 30 after medication, and
## the respective counts
##
aux <- CountBeforeAfterMedicationPerParticipant(x = scases, beforeThr = 30, afterThr = 30)
sel <- aux$sel
nba <- aux$nba

## organize the outputs of the personalized 
## classifier tests into a list
##
algos <- list(randomForest = o1, elasticNet = o4, logisticRegr = o3, extraTrees = o2)

## create a matrix with adjusted p-values (using
## Benjamini-Hochberg multiple testing correction),
## with rows indexing patients and columns indexing
## classifiers
##
apvals <- GetPvaluesAcrossAlgosAndParticipants(algos, "BH")$apvals

## adjust the UI-tests p-values using Benjamini-Hochberg
## correction, and catenate the adjusted p-values from 
## the classifier and UI tests
##
aui <- apply(o0[[1]], 2, p.adjust, "BH")
apvals <- cbind(apvals, aui)

## create a matrix with median AUROC values,
## with rows indexing patients and columns indexing
## classifiers
##
auroc <- GetAurocAcrossAlgosAndParticipants(algos)

## sort the participants according to the median AUROC 
## from the random forests classifier
##
oParticipants <- names(sort(o1[[1]][, 1], decreasing = TRUE))

## generate Figure 8
##
mat <- matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3), 7, 3, byrow = TRUE)
postscript(paste(fig.path, "fig8.ps", sep = ""), width = 10, height = 8)
layout(mat)
par(mar = c(0, 4, 0, 2) + 0.1, mgp = c(2.5, 1, 0)) #mar = c(bottom, left, top, right)
dat <- nba
dat <- dat[match(oParticipants, rownames(dat)),]
barplot(t(dat), beside = TRUE, las = 2, col = c("red", "blue"), xaxt = "n", 
        ylab = "sample size", cex.lab = 1.5)
legend("topright", legend = c("before", "after"), text.col = c("red", "blue"), 
       bty = "n", cex = 1.3)
text(-2.5, 150, "(a)", cex = 2)
AurocPlotAcrossParticipants(auroc, jit = 0.1,
                            legendPos = "topright", ylab = "AUROC", cex.lab = 1.5,
                            main = "", mymar = c(0, 4, 0, 2) + 0.1, cex = 2,
                            oParticipants = oParticipants,
                            cex.legend = 1.3)
text(-0.3, 0.9, "(b)", cex = 2)
PvaluePlotAcrossParticipants(apvals, jit = 0.1, cex.participant = 0.61, 
                             cex = 2, cex.lab = 1.5,
                             legendPos = "topright", ylab = "",
                             main = "", mymar = c(17, 4, 0, 2) + 0.1,
                             oParticipants = oParticipants,
                             cex.legend = 1.3)
text(-0.3, 29, "(c)", cex = 2)
abline(h = 3, col = "grey", lty = 1)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0))
dev.off()


########################################
## Figure 9
########################################

## load the outputs of the real data analysis
##
load("real_data_analysis_outputs.RData")
#load(synGet("syn4677553")@filePath) ## or load directly from Synapse

## sort the participants according to the median AUROC 
## from the random forests classifier, and select the
## top and bottom 3 patients
##
oParticipants <- names(sort(o1[[1]][, 1], decreasing = TRUE))
topIds <- oParticipants[1:3]
bottomIds <- oParticipants[55:57]

## check the top features for the top 3 patients
## (according to the random forest importance
## measure)
##
featnms <- colnames(scases)[9:32]
GetParticipantImp(o1$participantOutputs, id = topIds[1], featnms)
GetParticipantImp(o1$participantOutputs, id = topIds[2], featnms)
GetParticipantImp(o1$participantOutputs, id = topIds[3], featnms)

## check the top features for the top 3 patients
## (according to the individual t-test p-values)
##
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[1], type = c("ttest"))
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[2], type = c("ttest"))
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[3], type = c("ttest"))

## check the top features for the top 3 patients
## (according to the individual Wilcoxon rank
## sum test p-values)
##
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[1], type = c("wilcox"))
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[2], type = c("wilcox"))
GetUiTestTopFeatures(o0$participantOutputs, id = topIds[3], type = c("wilcox"))

## note that numberTaps, meanTappingInter, and medianTappingInter
## ranked among the top 3 features in all cases

## generate Figure 9
##
cl <- 1.0
cd <- 1.0
lcol <- "grey"
ylim1 <- c(20, 240)
ylim2 <- c(0.05, 0.7)
feat1 <- "numberTaps"
feat2 <- "meanTappingInter"
postscript(paste(fig.path, "fig9.ps", sep = ""), width = 10, height = 6)
par(mfrow = c(4, 3), mar = c(2, 2.5, 1, 0.5) + 0.1, mgp = c(1.5, 0.5, 0))
PlotBeforeAfterMedication(healthCodeId = topIds[1], featName = feat1, 
                          main = topIds[1], ylim = ylim1,
                          dat = scases, ylab = "number of taps", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(a)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = topIds[2], featName = feat1, 
                          main = topIds[2], ylim = ylim1, 
                          dat = scases, ylab = "number of taps", xlab = "", 
                          cex.dot = cd, lcol = lcol)
mtext("(b)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = topIds[3], featName = feat1, 
                          main = topIds[3], ylim = ylim1, 
                          dat = scases, ylab = "number of taps", xlab = "", 
                          cex.dot = cd, lcol = lcol)
mtext("(c)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = topIds[1], featName = feat2, 
                          main = topIds[1], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "", 
                          cex.dot = cd, lcol = lcol)
mtext("(d)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = topIds[2], featName = feat2, 
                          main = topIds[2], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "", 
                          cex.dot = cd, lcol = lcol)
mtext("(e)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = topIds[3], featName = feat2, 
                          main = topIds[3], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "", 
                          cex.dot = cd, lcol = lcol)
mtext("(f)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
###
PlotBeforeAfterMedication(healthCodeId = bottomIds[1], featName = feat1, 
                          main = bottomIds[1], ylim = ylim1, legendPos = "bottomright",
                          dat = scases, ylab = "number of taps", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(g)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = bottomIds[2], featName = feat1, 
                          main = bottomIds[2], ylim = ylim1, legendPos = "bottomright", 
                          dat = scases, ylab = "number of taps", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(h)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = bottomIds[3], featName = feat1, 
                          main = bottomIds[3], ylim = ylim1, legendPos = "bottomright", 
                          dat = scases, ylab = "number of taps", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(i)", side = 2, at = 230, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = bottomIds[1], featName = feat2, 
                          main = bottomIds[1], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(j)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = bottomIds[2], featName = feat2, 
                          main = bottomIds[2], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(k)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
PlotBeforeAfterMedication(healthCodeId = bottomIds[3], featName = feat2, 
                          main = bottomIds[3], ylim = ylim2, 
                          dat = scases, ylab = "mean tapping interval", xlab = "",
                          cex.dot = cd, lcol = lcol)
mtext("(l)", side = 2, at = 0.67, cex = cl, line = -2, las = 2)
par(mfrow = c(1, 1), mgp = c(3, 1, 0), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

