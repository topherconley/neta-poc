
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
  nodelist <- getnodelist(maxpernode=ncores, f = "nodelist_01.txt")
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
filterdir <- "/home/cconley/scratch-data/neta-poc/nohc-mad-filtered/"
fy <- "poc-mrna-eset-dropout-std.rds"
fx <- "poc-cna-eset-dropout-std.rds"
library(Biobase)
xset <- readRDS(file.path(filterdir, fx))
X <- t(exprs(xset))
yset <- readRDS(file.path(filterdir, fy))
Y <- t(exprs(yset))
#assert the sample names match
stopifnot(all(sampleNames(xset) == sampleNames(yset)))

##################
#PREVIOUS TUNING #
##################

#prev <- new.env()
#tunedir <- "~/scratch-data/neta-poc/mrna-tuning/05"
#load(file = file.path(tunedir, "poc-spacemap-05.rda"), envir = prev)
#tune <- prev$cvsmap$minTune
tune <- list(lam1 = 132, lam2 = 40, lam3 = 5) 

##################
#    ENSEMBLE    #
##################

#result directory
respath <- "/home/cconley/scratch-data/neta-poc/mrna-boot-vote/01"

library(spacemap)
seed <- 615621
tictoc <- system.time({ens <- spacemap::bootEnsemble(Y = Y, X = X, tune = tune,
                                                     method = "spacemap", B = 200,
                                                     resPath = respath,
                                                     seed = seed, p0 = 1,
                                                     tol = 1e-4, cdmax = 90e7)})
save.image(file = file.path(respath, "poc-boot-vote-p100-B200.rda"))
#stop cluster to not run out of memory
stopCluster(cl)
#re-register the sequential backend
registerDoSEQ()
object.size(ens)
bv <- bootVote(ens)

##################
#  SAVE RESULTS  #
##################
saveRDS(object = bv, file = file.path(respath, "poc-boot-vote-p100-B200.rds"))
save.image(file = file.path(respath, "poc-boot-vote-p100-B200.rda"))
