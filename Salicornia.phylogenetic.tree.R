library(data.table)
library(ape)
library(phyclust)
library(tidyverse)
library(ggplot2)
library(ggExtra)
library(gridExtra)
options(scipen=999)
library(data.table)


################################### quality filtered k-mer loci########################  ###################


geno3 <- read.delim("kmer_matrix.319.salicornia_accessions.filtered.5.percent.maf.txt", header = T, sep="\t")
dim(geno3)


rm_col <- c("SalD_143b")## remove SalD_143b, the duplicated.lines


geno3 <- geno3[, !(colnames(geno3) %in% rm_col), drop=FALSE] 
dim(geno3)


line_information = read.delim("line.info.salicornia2025.txt", header = TRUE, sep = "\t", quote = "\"") ##
dim(line_information)

geno_file = t(as.matrix( geno3[, 1:ncol(geno3)] ) ) # only genotypes
dim(geno_file)

## calculate genetic distances
ifelse(test = file.exists('./data1/distMat.RData'),
       yes = load('./data1/distMat.RData'), no = distMat <- dist(geno_file))
##  saving distance matrix
if (!file.exists('./data1/distMat.RData')) save(distMat, file = "./data1/distMat.RData")

load('./data1/distMat.RData')
hc2 <- as.phylo(hclust(distMat))
hc2$tip.label
hc2$edge

edgecols <- cbind('SalD'=NA, 1:nrow(hc2$edge), color='black') # create data frame
edgecols[, 1] <- hc2$tip.label[hc2$edge[, 2]] # 
edgecols <- as.matrix(merge(edgecols, line_information, by = 'SalD', all.x = T)) # get samples info
edgecols <- edgecols[order(as.numeric(edgecols[, 2])), ] # get samples in original order
edgecols[,4] <- as.numeric(edgecols[,4])

edgecols[, 3][edgecols[, 4] == "1"] =	"#6122b9ff" ## europaea_diploid israel 
edgecols[, 3][edgecols[, 4] == "2"] =	"#6172b7ff" ## sinus persica
edgecols[, 3][edgecols[, 4] == "3"] =	"#a9c8e1ff" ## persica iranica
edgecols[, 3][edgecols[, 4] == "4"] =	"#F7DADF" ## brachita
edgecols[, 3][edgecols[, 4] == "5"] =	"#8f4d69ff" ## europaea_tetraploid
edgecols[, 3][edgecols[, 4] == "6"] =	"lightsalmon1" ## bigelovii
edgecols[, 3][edgecols[, 4] == "7"] =	"khaki2" ##  S. species from korea
edgecols[, 3][edgecols[, 4] == "8"] =	"#FDAE61" ## veneta
edgecols[, 3][edgecols[, 4] == "9"] =	"#ECACB9" ## fructicosa
edgecols[, 3][edgecols[, 4] == "10"] =	"#BAE4BC" ##all.techticornia
edgecols[, 3][edgecols[, 4] == "18"] =  "gray30" ## batch7_korea
edgecols[, 3][edgecols[, 4] == "19"] =  "#C8A6B4" ## tetraploid S. procumbens
edgecols[, 3][edgecols[, 4] == "20"] =  "#6172b7ff" ## hybrid
edgecols[, 3][edgecols[, 4] == "21"] =  "brown2" ## batch10_europaea
edgecols[, 3][edgecols[, 4] == "22"]  = "#350571"   ##  europaea UAE 
edgecols[, 3][edgecols[, 4] == "23"]  = "magenta"   ##  europaea Egypt



tipcols <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))


  #plotting tree
pdf("output3/figure4.manuscript.pdf", width = 11.5, height = 7.5)
plotnj(unrooted.tree = hc2, type = 'u',
    show.tip.label = F, lab4ut = "axial", label.offset = 1, cex = 0.65,
    edge.color = edgecols[, 3], edge.width = 1, tip.color = tipcols[, 3], rotate.tree = 320)

dev.off()


## circular.phylogenetic tree

pdf("output3/figure.S9.manuscript.circular-phylogenetic.tree.pdf", width = 11.5, height = 7.5)
plotnj(unrooted.tree = hc2, type = 'fan',
       show.tip.label = T, lab4ut = "axial", label.offset = 1, cex = 0.38,
       edge.color = edgecols[, 3], edge.width = 1, tip.color = tipcols[, 3], rotate.tree = 320)
legend("topleft", lty=1, lwd = 2, cex = 0.65, box.lty =0,
   legend = c("europaea_israel", "sinus persica",	"persica iranica", "brachita", "europaea_tetraploid", "bigelovii", "herbaceae-like", "veneta", "fructicosa", "all.tecticornia", "batch7_korea", "S. procumbens", "hybrid", "batch10_europaea", "europaea_UAE", "europaea_egypt"),

    text.col = 'black', col = c("#6122b9ff", "#6172b7ff", "#a9c8e1ff", "#F7DADF", "#8f4d69ff", "lightsalmon1", "khaki2", "#FDAE61", "#ECACB9", "#BAE4BC",  "gray30", "#C8A6B4", "#6172b7ff", "brown2", "#350571", "magenta")) ## should be europaea israel, sinus perisa and persica iranica
dev.off()




# compute A matrix for eigen values and plot PCA
library(rrBLUP)
geno_file = t( as.matrix( geno3[, 1:ncol(geno3)] ) ) # only genotypes
geno_file[1:15,1:5]
A = A.mat(geno_file, impute.method="mean", return.imputed = T)

# PCA label names
names = as.data.frame(rownames(A$A)); colnames(names)=c("SalD")
names = merge(names, edgecols, by='SalD', sort = F, all.x = T)

##write.table(names, file = "Salicornia.names.updated1.SalD.kmer-based.txt", row.names = F, quote = F, sep = '\t')


# give shapes to points
names$color = as.character(names$color)
names$subspecies = as.numeric(names$subspecies)

# eigenvectors
library(rrBLUP)
e = eigen(A$A)

## PCA plot
pdf(file = 'output3/PCA.plot.Figure.S10.pdf', width = 11.25, height = 8)

#par(mfrow=c(2,1)) ## this is for two plot in a same page
plot(e$vectors[,1], e$vectors[,2], pch=names$subspecies, cex=1.05, col=names$color, lwd =1.5,
     xlab=paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''),
     ylab=paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''),
     main = substitute(paste("PCA of dicoccides species")),
     font.main=1, cex.main=1.5, cex.lab=0.75)
legend("topright", pch=c(1,2,3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21), inset=c(0, 0), cex=0.95, pt.cex = 1.85, pt.lwd = -2, box.lty =0,

col = c("#6122b9ff", "#6172b7ff", "#a9c8e1ff", "#F7DADF", "#8f4d69ff", "lightsalmon1", "khaki2", "#FDAE61", "#ECACB9", "#BAE4BC",  "gray30", "#C8A6B4", "#6172b7ff", "brown2", "#350571", "magenta" ),
       legend = c("europaea_israel", "sinus persica",	"persica iranica", "brachita", "europaea_tetraploid", "bigelovii", "herbaceae-like", "veneta", "fructicosa", "all.tecticornia", "batch7_korea", "S. procumbens", "hybrid", "batch10_europaea", "europaea_UAE", "europaea_egypt"))
dev.off()










