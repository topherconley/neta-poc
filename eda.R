

library(readr)
full <- read_csv(file = "~/scratch-data/neta-poc/cell-paper-table-s2-data/FullProteomeOvarianData.csv")
vpoc <- read_csv(file = "~/scratch-data/neta-poc/cell-paper-table-s2-data/VariableProteomeOvarianData.csv")
vpoc_samples <- colnames(vpoc)[3:length(vpoc)]

#Divide samples by collection site and reduce to 169 from 174
samples <- colnames(full)[3:length(full)]
tcga_samples <- gsub("JHU-|PNNL-", "",samples)
not_hgsc <- tcga_samples[!(tcga_samples %in% vpoc_samples)]
not_hgsc_index <- sapply(not_hgsc, function(pat) grep(pat = pat, x = samples))
samples <- samples[-not_hgsc_index]
pnnl_samples <- samples[grep("PNNL", samples)]
jhu_samples <- samples[grep("JHU", samples)]
tcga_jhu_samples <- gsub("JHU-", "", jhu_samples)
tcga_pnnl_samples <- gsub("PNNL-", "", pnnl_samples)
intersect(tcga_jhu_samples, tcga_pnnl_samples)

mjhu <- as.matrix(full[,jhu_samples])
rownames(mjhu) <- full$refseq_peptide
mpnnl <- as.matrix(full[,pnnl_samples])
rownames(mpnnl) <- full$refseq_peptide

########################
library(Biobase)
ed <- as.matrix(vpoc[,vpoc_samples])
rownames(ed) <- vpoc$refseq_peptide
#initial feature information
fd <-  data.frame( RefSeqProteinID = vpoc$refseq_peptide, 
                   row.names = vpoc$refseq_peptide)
vpocset <- ExpressionSet(assayData = ed,
                         featureData = AnnotatedDataFrame(data = fd))
pexp <- exprs(vpocset)
#no missing proteins
stopifnot(!any(is.na(pexp)))
#Gene-wise mean-sd 

sd_pexp <- apply(X = pexp, MARGIN = 1, FUN = sd)
mad_pexp <- apply(X = pexp, MARGIN = 1, FUN = mad)
mean_pexp <- apply(X = pexp, MARGIN = 1, FUN = mean)

#subset of full to match the variable set
mjhu_fullv <- mjhu[featureNames(vpocset),]
mpnnl_fullv <- mpnnl[featureNames(vpocset),]

lsd_sample <- list(variable = apply(X = pexp, MARGIN = 2, FUN = sd),
                   full_jhu = apply(X = mjhu, MARGIN = 2, FUN = sd, na.rm = T),
                   var_jhu = apply(X = mjhu_fullv, MARGIN = 2, FUN = sd, na.rm = T),
                   full_pnnl = apply(X = mpnnl, MARGIN = 2, FUN = sd, na.rm = T),
                   var_pnnl = apply(X = mpnnl_fullv, MARGIN = 2, FUN = sd, na.rm = T))
lsd_gene <- list(variable = apply(X = pexp, MARGIN = 1, FUN = sd),
                 full_jhu = apply(X = mjhu, MARGIN = 1, FUN = sd, na.rm = T),
                 var_jhu = apply(X = mjhu_fullv, MARGIN = 1, FUN = sd, na.rm = T),
                 full_pnnl = apply(X = mpnnl, MARGIN = 1, FUN = sd, na.rm = T),
                 var_pnnl = apply(X = mpnnl_fullv, MARGIN = 1, FUN = sd, na.rm = T))
lmean_sample <- list(variable = apply(X = pexp, MARGIN = 2, FUN = mean),
                   full_jhu = apply(X = mjhu, MARGIN = 2, FUN = mean, na.rm = T),
                   var_jhu = apply(X = mjhu_fullv, MARGIN = 2, FUN = mean, na.rm = T),
                   full_pnnl = apply(X = mpnnl, MARGIN = 2, FUN = mean, na.rm = T),
                   var_pnnl = apply(X = mpnnl_fullv, MARGIN = 2, FUN = mean, na.rm = T))
lmean_gene <- list(variable = apply(X = pexp, MARGIN = 1, FUN = mean),
                 full_jhu = apply(X = mjhu, MARGIN = 1, FUN = mean, na.rm = T),
                 var_jhu = apply(X = mjhu_fullv, MARGIN = 1, FUN = mean, na.rm = T),
                 full_pnnl = apply(X = mpnnl, MARGIN = 1, FUN = mean, na.rm = T),
                 var_pnnl = apply(X = mpnnl_fullv, MARGIN = 1, FUN = mean, na.rm = T))

gggene <- data.frame(mean = unlist(lmean_gene), sd = unlist(lsd_gene), 
                     source = factor(rep(names(lmean_gene), times = sapply(lmean_gene, length)), 
                                     levels = c("variable", "full_jhu", "var_jhu", "full_pnnl", "var_pnnl")), 
                     margin = rep("gene", length(unlist(lmean_gene))))


ggsample <- data.frame(mean = unlist(lmean_sample), sd = unlist(lsd_sample), 
                       source = factor(rep(names(lmean_sample), times = sapply(lmean_sample, length)), 
                                       levels = c("variable", "full_jhu", "var_jhu", "full_pnnl", "var_pnnl")), 
                       margin = rep("sample", length(unlist(lmean_sample))))


gggene_var <- gggene[grep("var", gggene$source),]
ggsample_var <- ggsample[grep("var", ggsample$source),]


#subset vs full comparision 
library(ggplot2)
ggplot(data = gggene, aes(x = mean, y = sd, colour = source)) + geom_point(alpha = 0.3)  + theme_bw() + 
  geom_smooth(se = T)
ggplot(data = ggsample, aes(x = mean, y = sd, colour = source)) + geom_point(alpha = 0.3)  + theme_bw() + 
  geom_smooth(se = T)
ggplot(data = gggene, aes(x = mean,colour= source)) + geom_density() + theme_bw() 

#subset vs full(subset) comparison
library(ggplot2)
ggplot(data = gggene_var, aes(x = mean, y = sd, colour = source)) + geom_point(alpha = 0.3)  + theme_bw() + 
  geom_smooth(se = T)
ggplot(data = ggsample_var, aes(x = mean, y = sd, colour = source)) + geom_point(alpha = 0.3)  + theme_bw() + 
  geom_smooth(se = T)
ggplot(data = gggene_var, aes(x = mean,colour= source)) + geom_density() + theme_bw() 




####missing values



#percent missing genes per sample
summary(apply(mjhu, 2, function(x) sum(is.na(x))/dim(mjhu)[1]))
summary(apply(mpnnl, 2, function(x) sum(is.na(x))/dim(mpnnl)[1]))


#percent missing samples per gene
summary(apply(mjhu, 1, function(x) sum(is.na(x))/dim(mjhu)[2]))
summary(apply(mpnnl, 1, function(x) sum(is.na(x))/dim(mpnnl)[2]))


sum(apply(mjhu, 1, function(x) all(!is.na(x))) & apply(mpnnl, 1, function(x) all(!is.na(x))))

mjhu2 <- mjhu[!(apply(mjhu, 1, function(x) sum(is.na(x))/dim(mjhu)[2]) > 0.20),]
dim(mjhu2)
dim(mjhu)

mpnnl2 <- mpnnl[!(apply(mpnnl, 1, function(x) sum(is.na(x))/dim(mpnnl)[2]) > 0.20),]
dim(mpnnl)
dim(mpnnl2)

#isoform issue

cmn_genes <- intersect(rownames(mjhu2), rownames(mpnnl2))
sum(rownames(mjhu2) %in% cmn_genes)
sum(rownames(mpnnl2) %in% cmn_genes)


#cmn_samples <- intersect(colnames(mpnnl3), colnames(mjhu3))
#mjhu3 <- mjhu2[cmn_genes,cmn_samples]
#mpnnl3 <- mpnnl2[cmn_genes,cmn_samples]

mjhu3 <- mjhu2[cmn_genes,]
mpnnl3 <- mpnnl2[cmn_genes,]


library(impute)
mjhu4 <- impute::impute.knn(data = mjhu3)$data
mpnnl4 <- impute::impute.knn(data = mpnnl3)$data

mcmn <- rbind(t(mpnnl4), t(mjhu4))
dim(mcmn)
#rownames(mcmn) <- paste0(rep(c("pnnl", "jhu"), each = 24), "-",rownames(mcmn))


#collection_site <- paste0(rep(c("PNNL", "JHU"), times  = c(84,122)), "-",rownames(mcmn))
#rownames(mcmn) <- collection_site
pca <- prcomp(x = mcmn, scale. = TRUE)
screeplot(pca)
#library(ggfortify)
#autoplot(pca, data = data.frame(site = rep(c("pnnl", "jhu"), each = 24)), colour = "site") + theme_bw()

#ggiris <- cbind(as.data.frame(pca$x), site = data.frame(site = rep(c("pnnl", "jhu"), each = 24)))
ggiris <- cbind(as.data.frame(pca$x), site = data.frame(site = rep(c("pnnl", "jhu"), times  = c(ncol(mpnnl4), ncol(mjhu4)))))
library(ggplot2)
gg12 <- ggplot(data = ggiris, aes(x=PC1, y=PC2, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()
gg13 <- ggplot(data = ggiris, aes(x=PC1, y=PC3, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()
gg23 <- ggplot(data = ggiris, aes(x=PC2, y=PC3, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()
gg14 <- ggplot(data = ggiris, aes(x=PC1, y=PC4, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()
gg24 <- ggplot(data = ggiris, aes(x=PC2, y=PC4, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()
gg34 <- ggplot(data = ggiris, aes(x=PC3, y=PC4, group = site, colour = site)) + geom_point(size = 2) + scale_color_brewer(type = "qual") + theme_bw()

library(gridExtra)
grid.arrange(grobs = list(gg12, gg13, gg23, gg14, gg24, gg34), ncol = 2)





