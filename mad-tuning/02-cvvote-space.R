####################
# PARALLEL BACKEND #
####################


#Take in arguments from sbatch. 
args <- commandArgs(TRUE)

if (length(args) == 0) { 
  ncores <- 1
  comp <- "singlenode"
} else {
  ncores <- as.integer(args[1])
  comp <- "gauss" 
}

if (comp == "gauss") {
  
  #====Setup for running on Gauss... ======#
  
  getnodelist <- function(maxpernode=30, f = "nodelist.txt") {
    nodelist <- Sys.getenv('SLURM_NODELIST')
    system(paste0("scontrol show hostname ", nodelist, " > ", f))
    x <- if (nzchar(nodelist)) readLines(f) else rep('localhost', 3)
    d <- as.data.frame(table(x), stringsAsFactors=FALSE)
    rep(d$x, pmax(d$Freq, maxpernode))
  }
  
  library(doParallel)
  nodelist <- getnodelist(maxpernode=ncores, f = "nodelist_02.txt")
  print(nodelist)
  cl <- makePSOCKcluster(nodelist, outfile='')
  registerDoParallel(cl)
} else if (comp == "singlenode") {
  #===== Setup for running on a multicore machine =====#
  suppressPackageStartupMessages(library(doParallel))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  stop("Argument: comp is not specified correctly.")
}


##################
#     Data       #
##################
filterdir <- "/home/cconley/scratch-data/neta-poc/mad-filtered-s2-data"
#filterdir <- "~/scratch-data/neta-poc/mad-filtered-s2-data"
fy <- "poc-prot-eset-dropout-std.rds"
fx <- "poc-cna-eset-dropout-std.rds"
library(Biobase)
xset <- readRDS(file.path(filterdir, fx))
X <- t(exprs(xset))
yset <- readRDS(file.path(filterdir, fy))
Y <- t(exprs(yset))
#assert the sample names match
stopifnot(all(sampleNames(xset) == sampleNames(yset)))

##################
#   INIT FIT     #
##################
N <- nrow(Y)
Q <- ncol(Y)
lam1start <- function(n, q, alpha) { 
  sqrt(n) * qnorm(1 - (alpha/ (2*q^2)))
}
lam0 <- lam1start(n = floor(N - N*.10), q = Q, alpha = 1e-5)
lam0

##################
#   TUNING       #
##################

#training/test setse
trainSets <- readRDS(file.path(filterdir, "train_sets.rds"))
testSets <- readRDS(file.path(filterdir, "test_sets.rds"))

##################
#   TRY O1       #
##################

#GRID
tmap <- expand.grid(lam1 = exp(seq(log(150), log(250), length = 100)))

#result directory
respath <- "/home/cconley/scratch-data/neta-poc/mad-tuning/02"
if (!dir.exists(respath)) { 
  system(paste("mkdir -p", respath))  
}

library(spacemap)
tictoc <- system.time({cvsmap <- spacemap::cvVote(Y = Y,
                                                  trainIds = trainSets, testIds = testSets, 
                                                  method = "space", tuneGrid = tmap, 
                                                  resPath = respath,
                                                  tol = 1e-4, cdmax = 90e7)})
save.image(file = file.path(respath, "poc-space-02.rda"))
stopCluster(cl)
