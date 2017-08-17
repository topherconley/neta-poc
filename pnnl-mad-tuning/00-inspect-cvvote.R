outdir <- "~/scratch-data/neta-poc/pnnl-mad-tuning/01"
load(file = file.path(outdir, "poc-space-01.rda"))

outdir <- "~/scratch-data/neta-poc/pnnl-mad-tuning/02"
load(file = file.path(outdir, "poc-space-02.rda"))

outdir <- "~/scratch-data/neta-poc/pnnl-mad-tuning/03"
load(file = file.path(outdir, "poc-space-03.rda"))


library(ggplot2)
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

cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")
library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)


library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)

cvsmap$logcvScore
dim(cvsmap$cvVote$yy)


summary(tmap[cvsmap$metricScores$dfGamma[,1] > 500 & cvsmap$metricScores$dfGamma[,1] < 10000,"lam2"])


tmap <- expand.grid(lam1 = c(182.141,seq(150, 200, length = 5)),
                    lam2 = sqrt(seq(10^2, 70^2, length = 10)),
                    lam3 = (seq(sqrt(0), sqrt(90), length = 10))^2)


tmap <- expand.grid(lam1 = seq(sqrt(120), sqrt(200), length = 10)^2,
                    lam2 = seq(30, 60, length = 10),
                    lam3 = seq(0, 100, length = 10))


