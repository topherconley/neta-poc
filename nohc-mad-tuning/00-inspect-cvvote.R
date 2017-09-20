outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/01"
load(file = file.path(outdir, "poc-space-01.rda"))

outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/02"
load(file = file.path(outdir, "poc-spacemap-02.rda"))

outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/03"
load(file = file.path(outdir, "poc-spacemap-03.rda"))

outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/04"
list.files(path = outdir, "drop")
load(file = file.path(outdir, "poc-spacemap-04-dropNA.rda"))


outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/05"
load(file = file.path(outdir, "poc-spacemap-05.rda"))


outdir <- "~/scratch-data/neta-poc/nohc-mad-boot-vote/02/"
load(file = file.path(outdir, "poc-boot-vote-p90-B1000.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p90-B1000.rds"))

outdir <- "~/scratch-data/neta-poc/nohc-mad-boot-vote/03/"
#load(file = file.path(outdir, "poc-spacemap-05.rda"))
load(file = file.path(outdir, "poc-boot-vote-p100-B1000.rda"))
bv <- readRDS(file = file.path(outdir, "poc-boot-vote-p100-B1000.rds"))



library(ggplot2)
outdir <- "~/scratch-data/neta-poc/nohc-mad-tuning/03"
load(file = file.path(outdir, "poc-spacemap-03.rda"))

nsplits <- sapply(testSets, length)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)


cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)

cvsmap$metricScores$rss[cvsmap$minIndex,]
mean(cvsmap$metricScores$dfParCor[cvsmap$minIndex,])
mean(cvsmap$metricScores$dfGamma[cvsmap$minIndex,])

cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")
library(gridExtra)
grid.arrange(cvVis1[[1]] + coord_cartesian(ylim = c(7.105, 7.111)), cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)



#########################
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
cvVis3 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")

grid.arrange(cvVis1[[1]] + coord_cartesian(ylim = c(7.106, 7.12)),
             cvVis1[[2]],
             cvVis2[[1]] + coord_cartesian(ylim = c(7.106, 7.12)),
             cvVis3[[1]] + coord_cartesian(ylim = c(7.106, 7.12)),
             ncol = 2)


grid.arrange(cvVis1[[1]], cvVis2[[1]],cvVis3[[1]], ncol = 2)

tmap <- expand.grid(lam1 = c(131.7571, seq(131, 132.2, length = 5)),
                    lam2 = c(44.50996, seq(43.9, 45, length = 5)),
                    lam3 = seq(4, 8, length = 6))

library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)

library(spacemap)
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)




library(ggplot2)
degxy <- rowSums(cvsmap$cvVote$xy)
degxy <- rowSums(bv$bv$xy)
qplot(degxy[degxy != 0]) + geom_histogram() + theme_bw()

library(ggplot2)
cvyy <- cvsmap$cvVote$yy
cvyy <- bv$bv$yy
cvyy[upper.tri(cvyy, diag = T)] <- 0
degyy <- rowSums(cvyy)
qplot(degyy[degyy != 0]) + geom_histogram() + theme_bw()
qplot(degxy[degxy != 0]) + geom_histogram() + theme_bw()


sum(rowSums(cvsmap$cvVote$xy) > 8)
quantile(rowSums(cvsmap$cvVote$xy), probs = seq(0.5, 1, length = 10))

cvsmap$logcvScore
dim(cvsmap$cvVote$yy)
cvsmap$minTune


summary(tmap[cvsmap$metricScores$dfGamma[,1] > 500 & cvsmap$metricScores$dfGamma[,1] < 10000,"lam2"])



