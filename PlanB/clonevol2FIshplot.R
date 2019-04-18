# dependencies of clonevol
library(gridBase)
library(gridExtra)
library(ggplot2)
library(igraph)
library(packcircles)
library(trees)

# dependencies of fishplot
library(png)
library(Hmisc)
library(plotrix)

library(clonevol)
library(fishplot)


### need to be automated



sc2 = read.table("/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/CCF_data2.txt",sep="\t",stringsAsFactors=F,header=T)
vafs2 = data.frame(sc2[,2]+1,sc2[,3:9]*100)
samples = c("t1", "t2", "t3", "t4", "t5","t6", "t7")
names(vafs2)[1] = "cluster"
names(vafs2)[2:8] = samples
vafs3 <- vafs2[which(vafs2$cluster %in% cluster_ls2[,1]),]

cluster_ls3 <- as.data.frame(table(vafs2[,1]))
cluster_ls2 <- cluster_ls3[which(cluster_ls3$Freq > 5), ]
vafs3[vafs3$cluster == 22, 1] <- 2
vafs3[vafs3$cluster == 106, 1] <- 3


#------(should be replaced by EvolCancer)------#
## run clonevol 
res = infer.clonal.models(variants=vafs3, cluster.col.name="cluster", ccf.col.names=samples,
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          cluster.center="mean", num.boots=1000, founding.cluster=1,
                          min.cluster.vaf=0.01, p.value.cutoff=0.01, alpha=0.1, random.seed=63108)

## create a list of fish objects - one for each model (in this case, there's only one)
f = generateFishplotInputs(results=res)
fishes = createFishPlotObjects(f)

## plot each of these with fishplot
pdf('fish.pdf', width=8, height=4)
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="PatientID", cex.title=0.7,
           vlines=seq(1, length(samples)), vlab=samples, pad.left=0.5)
}
dev <- dev.off()

