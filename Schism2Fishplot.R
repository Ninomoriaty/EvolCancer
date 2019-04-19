# dependencies of fishplot
library(png)
library(Hmisc)
library(plotrix)
library(fishplot)

## prepare Schism input
dir.cluster.tsv = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/hu.cluster.tsv"
dir.loci.tsv = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/hu.loci.tsv"
dir.output = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism"

prepareSchismInput <- function(dir.cluster.tsv, dir.loci.tsv, dir.output){
  # test: separate the first 3 sample and the next 4 sample to genrate a more proper result
  sep3 <- c("WGC033386D", "WGC033387D", "WGC033388D")
  sep4 <- c("WGC033389D", "WGC033391D", "WGC033392D", "WGC033393D")
  # read mutations in cluster.tsv from PyClone and get targeted clusters
  cluster.tsv = read.table(dir.cluster.tsv, sep="\t", stringsAsFactors=F, header=T)
  cluster_filter = cluster.tsv[!is.na(cluster.tsv$cluster_id) & cluster.tsv$size >= 5 & cluster.tsv$mean >= 0.1,]
  cluster_ls  = unique(cluster_filter$cluster_id)
  
  # read mutations within targeted clusters in loc.tsv file from PyClone 
  loci.tsv = read.table(dir.loci.tsv, sep="\t", stringsAsFactors=F, header=T)
  loci_filtered <- loci.tsv[which(loci.tsv$cluster_id %in% cluster_ls & loci.tsv$sample_id %in% sep4),]
  
  # generate two outputs for SCHISM: clusterEstimates.tsv, mutation_to_cluster.tsv
  clusterEstimates.tsv <- data.frame(loci_filtered$sample_id, loci_filtered$mutation_id, 
                                     loci_filtered$cellular_prevalence, loci_filtered$cellular_prevalence_std)
  names(clusterEstimates.tsv) = c("sampleID",	"mutationID",	"cellularity", "sd")
  
  mutation_to_cluster.tsv <- data.frame(loci_filtered$mutation_id, loci_filtered$cluster_id)
  names(mutation_to_cluster.tsv) = c("mutationID",	"clusterID")
  
  # set the directory of outputs and make two SCHISM-needed files 
  setwd(dir.output)
  write.table(clusterEstimates.tsv, file="W.clusterEstimates.tsv", sep = "\t", quote=F, row.names=F)
  write.table(mutation_to_cluster.tsv, file="W.mutation-to-cluster.tsv", sep = "\t", quote=F, row.names=F)
}

###### generate results for createFishPlotObjects ######
## Method2: get the result from GA.consensusTree, but it take so much time for runSchism to calculate the result.
dir.GA.consensusTree = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism/0_W_results.GA.consensusTree"
dir.cluster.cellularity = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism/0_W_results.cluster.cellularity"

schism2Fishplot <- function(dir.cluster.cellularity, dir.GA.consensusTree){
  # get cellularity infomation
  cluster.cellularity = read.table(dir.cluster.cellularity, sep="\t", stringsAsFactors=F, header=T)
  
  ls.sample = unique(cluster.cellularity$sampleID)
  ls.cluster = unique(cluster.cellularity$clusterID)
  
  frac.c = c()
  for (sample_name in ls.sample) {
    sample <- cluster.cellularity[which(cluster.cellularity$sampleID == sample_name), ]$cellularity
    frac.c <- c(frac.c,sample)
  }
  frac.table = matrix(frac.c*100, ncol = length(ls.sample))
  rownames(frac.table) <- ls.cluster
  colnames(frac.table) <- ls.sample
  
  
  ## generate tree information
  GA.consensusTree = read.table(dir.GA.consensusTree, sep="\t", stringsAsFactors=F, header=T)
  
  ls.ends = GA.consensusTree[which(!GA.consensusTree$parent %in% GA.consensusTree$child), 1] 
  
  # make sure the first node of the tree would be first in parents
  if (!rownames(frac.table)[1] %in% ls.ends){
    frac.end = frac.table[which(rownames(frac.table) == ls.ends[1]),]
    frac.table <- frac.table[-which(rownames(frac.table) == ls.ends[1]),]
    frac.table <- rbind2(frac.end, frac.table)
    rownames(frac.table)[1] <- ls.ends[1]
  }
  
  # figure out the subclonal evolution relationship
  parents <- c(1:length(ls.cluster))
  for (cluster_name in ls.cluster){
    if (cluster_name %in% ls.ends){
      parents[which(rownames(frac.table) == cluster_name)] <- 0
    } else{
      parents[which(rownames(frac.table) == cluster_name)] <- which(rownames(frac.table) == GA.consensusTree[which(GA.consensusTree$child == cluster_name), 1])
    }
  }
  
  # test: what if the value is ok 
  frac.table[3,1] <- 30
  frac.table[2,4] <- 10
  frac.table[3,5] <- 49.60
  frac.table[2,5] <- 49.50
  frac.table[2,6] <- 49.70
  
  # test: first 3 columns errors(couldn't make a combination of subclones)
  frac.table[3,1] <- 38.75 # minus 1
  
  # test: last 4 columns errors
  frac.table[4,2] <- 49.72
  frac.table[3,4] <- 43.00
  
  # tips: make sure that the parents start from 0 or you will get c stack error
  fish = createFishObject(frac.table, parents, timepoints = c(1:length(ls.sample)))
  fish = layoutClones(fish)
  fishPlot(fish, shape="spline", title.btm="PatientID", cex.title=0.5,
           vlines=seq(1, length(ls.sample)), vlab=ls.sample, pad.left=0.5)
}

### Method1: get the result from HT, but however, it's so strange that I couldn't find the relationship by cost.
## you need to know the method resolve this data, however, it is given by Schism..
# HT.cpov = read.table("/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism/0_W_results.HT.cpov", sep="\t", stringsAsFactors=F, header=T)
# parents_filter <- HT.cpov[which(HT.cpov$topologyCost == 0 & HT.cpov$parent.cluster != HT.cpov$child.cluster), c("parent.cluster", "child.cluster")]
# parents_freq <- as.data.frame(table(parents_filter$parent.cluster))
# order(as.data.frame(table(parents_filter$parent.cluster)), decreasing = T)
# parents_filter[order(, decreasing = T), ]

# 
# df.parents <- data.frame(parent.cluster = c(NA), child.cluster = c(NA))
# for (cluster_name in ls.cluster){
#   child_cluster <- parents_filter[which(parents_filter$parent.cluster == cluster_name), ]
#   if (length(child_cluster[,1]) == 0){
#     next()
#   } else {
#     df.parents <- df.parents[-which(df.parents$child.cluster %in% child_cluster$child.cluster),]
#   }
#   df.parents <- rbind(df.parents, child_cluster)
# }

## get parents
# ls.ends <- df.parents[-which(df.parents$parent.cluster %in% df.parents$child.cluster), ]$parent.cluster
# parents <- c(1:length(ls.cluster))
# 
# for (cluster_name in ls.cluster){
#   if (cluster_name %in% ls.ends){
#     parents[which(rownames(frac.table) == cluster_name)] <- 0
#   } else{
#     parents[which(rownames(frac.table) == cluster_name)] <- which(rownames(frac.table) == df.parents[which(df.parents$child.cluster == cluster_name), 1])
#   }
# }

