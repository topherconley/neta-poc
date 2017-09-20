
#SPACE
outdir <- "~/scratch-data/neta-poc/mrna-tuning/01"
load(file = file.path(outdir, "poc-space-01.rda"))
library(ggplot2)
nsplits <- sapply(testSets, length)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
outdir <- "~/scratch-data/neta-poc/mrna-tuning/02"
load(file = file.path(outdir, "poc-space-02.rda"))
cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
library(gridExtra)
grid.arrange(cvVis1[[1]], 
             cvVis2[[1]] ,ncol = 2)


#SPACEMAP
outdir <- "~/scratch-data/neta-poc/mrna-tuning/03"
load(file = file.path(outdir, "poc-spacemap-03.rda"))

outdir <- "~/scratch-data/neta-poc/mrna-tuning/04"
load(file = file.path(outdir, "poc-spacemap-04.rda"))
nsplits <- sapply(testSets, length)
#selecting high y--y, low x-y 
cvsmap$minTune
library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)

#visual inspection
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)

cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)


#find decent tuning parameter set that matches protein size
cvOut <- cvsmap
library(foreach)
avgMetrics <- function(cvOut, testSetLen) { 
  requireNamespace("ggplot2")
  total <- sum(testSetLen)
  foldconvmat <- !is.na(cvOut$metricScores[["rss"]])
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  requireNamespace("foreach")
  avgScores <- foreach(i = seq_along(cvOut$metricScores), .combine = 'cbind') %do% {
    score <- cvOut$metricScores[[i]]
    stype <- names(cvOut$metricScores)[[i]]
    msum <- apply(X = score , MARGIN = 1, FUN = sum, na.rm = TRUE) 
    
    ###Selection of tuning grid restricted to case when majority of test hold outs
    ###were evaluated (chunked by folds) 
    
    # divide rss by number of test sets with corresponding trained convergence
    #divide other metrics by the number of converged folds
    
    if (stype == "rss") { 
      ifelse(ntestconv > 0.5*total, msum / ntestconv, Inf)
    } else {
      ifelse(ntestconv > 0.5*total, msum / nfoldconv, Inf)
    }
  }
  colnames(avgScores) <- names(cvOut$metricScores)
  avgScores
}
mets <- avgMetrics(cvOut = cvsmap, testSetLen = sapply(testSets, length))
mets[,"rss"] <- log(mets[,"rss"])
tvis <- cbind(mets, tmap)

#looking at cv.vote drop
test_idx <- which(tvis$lam1 == 132 & tvis$lam2 == 40  & tvis$lam3 == 5)
tvis[test_idx,]
top_files <- list.files(path = "~/scratch-data/neta-poc/mrna-tuning/04/", pattern = "tuneid_040_", full.names = T)
source("~/repos/spacemap/R/cvVoteRDS.R")
top_vote <- cvVoteRDS(files = top_files, 
                      tol = 0.0,thresh = 0.5, method = "spacemap")
nonZeroUpper(top_vote$yy,0.0)
nonZeroWhole(top_vote$xy,0.0)
### 01 try 
#01 to get a ball park of boot.vote drop used 
#tvis[tvis$lam1 == 132 & tvis$lam2 == 40  & tvis$lam3 == 5,]
outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/01/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B200.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B200.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)

#################
#02 try
tvis <- tvis[tvis$lam1 > 134 & tvis$dfGamma > 2000,]
top_index <- tvis[,"rss"] < quantile(tvis[,"rss"], prob = 0.05) 
tvis_top <- tvis[top_index,]
tvis_top
#try this one
tvis_top[4,]

outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/02/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B200.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B200.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)

###############
#03 try
tvis <- tvis[tvis$lam1 > 134 & tvis$dfGamma > 3000,]
top_index <- tvis[,"rss"] < quantile(tvis[,"rss"], prob = 0.05) 
tvis_top <- tvis[top_index,]
tvis_top
#try this one
tvis_top[3,]

outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/03/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B200.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B200.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)

###############
#04 try
tvis <- tvis[tvis$lam1 > 134 & tvis$dfGamma > 5000,]
top_index <- tvis[,"rss"] < quantile(tvis[,"rss"], prob = 0.05) 
tvis_top <- tvis[top_index,]
tvis_top
#try this one
tvis_top[2,]

tvis <- tvis[tvis$dfGamma > 5000,]

outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/04/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B200.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B200.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)


#where are these in the global percentile? 
tvis <- cbind(mets, tmap)
quantile(x = tvis$rss, probs = seq(0, 1, by = 0.05))

#################################################
#SPACEMAP 05
outdir <- "~/scratch-data/neta-poc/mrna-tuning/05"
load(file = file.path(outdir, "poc-spacemap-05.rda"))
nsplits <- sapply(testSets, length)
#selecting high y--y, low x-y 
cvsmap$minTune
library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)

#visual inspection
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)

cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)


mets2<- avgMetrics(cvOut = cvsmap, testSetLen = sapply(testSets, length))
mets2[,"rss"] <- log(mets2[,"rss"])
tvis2 <- cbind(mets2, tmap)

###############
#04 try
tvis3 <- rbind(tvis, tvis2)
tvis3 <- tvis3[tvis3$lam1 == 135 & tvis3$dfGamma > 6000,]
top_index <- tvis3[,"rss"] < quantile(tvis3[,"rss"], prob = 0.05) 
tvis_top <- tvis3[top_index,]
tvis_top
#try this one
tvis_top[which.min(tvis_top$rss),]


###########################
###########################

#Final 
#cv vote drop
which(tvis2$lam2 == 32.6 & as.integer(tvis2$lam3) == 62)
tvis2[110,]
top_files <- list.files(path = "~/scratch-data/neta-poc/mrna-tuning/05/", pattern = "tuneid_110", full.names = T)
source("~/repos/spacemap/R/cvVoteRDS.R")
top_vote <- cvVoteRDS(files = top_files, 
                      tol = 0.0,thresh = 0.5, method = "spacemap")
# 2505 Y--Y
nonZeroUpper(top_vote$yy,0.0)
# 3624 X--Y 
nonZeroWhole(top_vote$xy,0.0)

outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/05/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B200.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B200.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)

outdir <- "~/scratch-data/neta-poc/mrna-boot-vote/06/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B1000.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B1000.rds"))
sum(sapply(ens, function(x) x$convergence))
library(spacemap)

# 1735 Y--Y 
nonZeroUpper(bv$bv$yy,0.0)
# 1231 X--Y 
nonZeroWhole(bv$bv$xy,0.0)

