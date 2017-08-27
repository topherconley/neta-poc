

library(devtools)
setwd("~/repos/spacemap/")
load_all()

mapNewIndex <- function(ccidx, minIndex) { 
  #map from na to cc indices
  downShift <- 0 
  newidx <- vector(mode = "integer", length = sum(ccidx))
  for(i in seq_along(ccidx)) { 
    if (!ccidx[i]) { 
      downShift <- downShift + 1
      next;
    }
    newidx[i - downShift] <- i - downShift
  }
  upi <- data.frame(old = which(ccidx), clean = newidx)
  upi$old[match(minIndex, upi$clean)]
}


dropNATune <- function(file) { 
  load(file)
  library(spacemap)
  ccidx <- complete.cases(cvsmap$metricScores$rss)
  if(any(!ccidx)) { 
    
    #skip if selected minIndex does not contain missing value
    if(!(cvsmap$minIndex %in% which(!ccidx))) { 
      return("Nothing dropped because NA not in CV-selected tuning folds.")
    }
    
    #remove tuning sets with NA
    metricScores <- lapply(cvsmap$metricScores, function(x) x[ccidx,])
    tmap <- tmap[ccidx,]
    #average across folds
    metricScoresAvg <- averageScores(metricScores = metricScores, testIds = testSets)
    
    ##Find the minimizing tuning parameter set
    newMinIndex <- minScoreIndex(metricScoresAvg)
    oldMinIndex <- mapNewIndex(ccidx = ccidx, minIndex = newMinIndex)
    
    if (metricScoresAvg[newMinIndex[1],"rss"] == Inf) { 
      stop("No convergence for any fold or any tuning parameter selection.")
    }   
    #obtain model fit files from best tune index
    tmp <- sub(pat = "/home/cconley", repl = "~/", x = respath)
    files <- list.files(path = tmp, 
                        pattern =  paste0("tuneid_", sprintf("%03d", oldMinIndex), "_fold"),
                        full.names = TRUE)
    voteFit <- cvVoteRDS(files = files, tol = 0.0, 
                         thresh = 0.5, method = "spacemap",
                         givenX = FALSE)
    minTune <- tmap[newMinIndex,]
    
    cvsmap <- list(cvVote = voteFit, minTune = minTune, minIndex = newMinIndex, 
                   logcvScore = log(metricScoresAvg[newMinIndex,"rss"]),
                   metricScores = metricScores, bestFoldFiles = files)
    dnaf <- paste0(strsplit(x = file, split = ".", fixed = T)[[1]][1],
                   "-dropNA.rda")
    save(list = ls(all.names = TRUE), file = dnaf, envir = environment())
  }
  
}

#
#
dropNATune("~/scratch-data/neta-poc/nohc-mad-tuning/04/poc-spacemap-04.rda")
