## First Attempt at analyzing SNP (vcf file) data

setwd("~/Documents/Postdoc/Adenocystis_biogeography/DataAnalyses/Adenocystis_Biogeography")

library(vcfR)
library(adegenet)
library(ape)

vcf <- read.vcfR("populations.snps.FirstAttempt.vcf")
head(vcf)
x <- vcfR2genlight(vcf)
x

x <- vcfR2genind(vcf)
x
sum(is.na(x$tab))

# create a PCA plot
X <- tab(x, freq = TRUE, NA.method = "mean")
class(X)

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

s.label(pca1$li)
title("PCA of Adeno dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)


s.class(pca1$li, pop(X))

# trying a chrom object?


#Extract a matric of Depth Coverage for each sample and each locus
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
head(dp)
dim(dp)
#Save coverage depth per individual per locus to file
write.csv(dp, file = "dp.csv")

#Calculate Mean Depth coverage for each locus
RowM <- rowMeans(dp[, -ncol(dp)],na.rm = TRUE)
RowM
mean(RowM)
range(RowM)
plot(RowM)


#Make histograms and Boxplots of Depth Coverage measurements
par(mfrow=c(2,2)) 
hist(RowM, breaks=100, xlab ="Mean Coverage of Each Loci", main = "Histogram of Mean Coverage \n of Each Loci")
boxplot(RowM,main ="Mean Coverage of Loci")


#Print the outliers
#Boxplot Outliers: outside 1.5 times the interquartile range above the upper quartile and below the lower quartile
outliersM = boxplot(RowM, plot=FALSE)$out
sort(outliersM)


### PCA plots

###LOAD LIBRARIES
library(adegenet)
library(vegan)
library(car)
library(rgl)
library(dartR)

my.genind <- vcfR2genind(vcf)
my.genind@all.names
head(indNames(my.genind),93)
locations<-read.csv("Locations.csv", header=FALSE)
locations
td <- read.table("Locations.csv", sep=',', colClasses="character")
td

pop(my.genind)<-factor(locations)
strata(my.genind) <- locations
my.genind
my.genind@strata
setPop(my.genind) <- ~V1
my.genind@pop

library(mmod)
diff_bov <- diff_stats(my.genind)
diff_bov$global

glPlot(x)



vcf <- read.vcfR('final.recode.vcf')
head(vcf)
head(getFIX(vcf))

#Extract the depth info
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
dp[1:4,1:6]

library(reshape2)
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)

palette=rep_len(c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87"), 93)

library(tidyverse)

ggplot(dpf, aes(x=Sample, y=Depth)) + geom_boxplot(fill=palette) + theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(trans=scales::log2_trans(), expand = c(0,0), breaks=c(1, 10, 100, 1000, 5000),
                     minor_breaks=c(1:10, 2:10*10, 2:10*100, 2:5*1000)) +
  theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6)) +
  theme(panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2))

##
heatmap.bp(dp[(sample(nrow(dp),1000)),], rlabels = FALSE, rbarplot = FALSE)
heatmap.bp(dp)
##


#Code for looking at missingness per sample
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)
myMiss <- data.frame((levels(dpf$Sample)), myMiss)
colnames(myMiss) <- c('Sample', 'Missing')
palette=rep_len(c("#F46D43", "#FDAE61", "#FEE08B", "#D9EF8B", "#A6D96A", "#66BD63"), 103)
ggplot(myMiss, aes(x=Sample, y=Missing)) + geom_col(fill=palette) + theme_bw() +
  labs(y='Missingness (%)') +theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=60,hjust=1))+scale_y_continuous(expand = c(0,0))



#make a NJ tree
tre <- njs(dist(as.matrix(GBS)))
plot(tre, "p", cex=0.5, FALSE, no.margin = TRUE, font=4, node.pos=1, edge.width=1.2)

#This gives a much nicer tree than with the iPyrad stuff

#To analyze PCA; first remove non-informative things
toRemove <- is.na(glMean(GBS, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
GBS2 <- GBS[, !toRemove]

glPlot(GBS2)

#Try reanalyzing with a NJ tree
tree2 <- njs(dist(as.matrix(GBS2)))
plot(tree2, "p", cex=0.5, FALSE, no.margin = TRUE, font=4, node.pos=1, edge.width=1.2)

#PCA plot
GBS_pca <- glPca(GBS, useC = FALSE, parallel = TRUE, nf=10)

#makes a bar plot of the eigenvalues and such
par(mar = c(4, 4, 4, 4))
barplot(GBS_pca$eig, main="eigenvalues", col=heat.colors(length(GBS_pca$eig)))

#Looks like the first 1 is very important; and the rest are not useful

#Now a high number of clusters are kept at first to find clusters
grp <- find.clusters(GBS, max.n.clust=40, glPca = GBS_pca, n.pca=7, choose = FALSE, stat = "BIC")

#and plot that
plot(grp$Kstat, type = "o", xlab = "number of clusters (K)") #again, as above, I normally set choose=TRUE and the plot is output automatically
     
#and plot that
plot(grp$Kstat, type = "o", xlab = "number of clusters (K)", ylab = "BIC", col = "blue",main = "Value of BIC versus number of clusters")
#again, as above, I normally set choose=TRUE and the plot is output automatically
     #So the BIC tells me that (1) it's way too zoomed in for some god awful reason and (2) that the 13th node is the lowest without including a lot of other junk in this dataset sooo we're going with that.
     
indNames(GBS2)

#Do a DAPC plot; n.pca keeps 100 levels because why not? n.da keeps 13 discriminate things
dapc_pops <- dapc(GBS2, glPca = GBS_pca, n.da = 10, n.pca = 6)
#Graph DAPC plot
scatter(dapc_pops, scree.pca = TRUE, bg="white", pch=20, cstar=0, col=palette, solid=.6,
        cex=3, clab=0, leg=TRUE, posi.pca="topleft")

#Reset the group for realz this time
grp <- find.clusters(GBS2, max.n.clust=12, glPca = GBS_pca, n.pca=100, n.clust = 12)

#Check to see what individuals correspond to what group
table(pop(GBS2), grp$grp)

#Reset the group for realz this time
grp <- find.clusters(GBS2, max.n.clust=3, glPca = GBS_pca, n.pca=100, n.clust = 3)

#Check to see what individuals correspond to what group
table(pop(GBS2), grp$grp)

vcf <- read.vcfR('populations.snps.vcf')
#vcf <- read.vcfR('populations.snps.FirstAttempt.vcf')
head(vcf)

GBS<- vcfR2genlight(vcf)

GBS@ind.names

factorname <- factor(c("BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North", "BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","Campbell Point","Campbell Point", "Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Watsons Road North","Watsons Road North","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South"))

factorname <- factor(c("BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","Campbell Point","Campbell Point", "Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Watsons Road North","Watsons Road North","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South"))

factorname <- factor(c("BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","Campbell Point","Campbell Point","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Watsons Road North","Watsons Road North","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South"))

length(factorname)

GBS@pop <- factorname
GBS@pop
poplevels<-GBS@pop

myFreq <- glMean(GBS)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="#5B88C1", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20,ylim = c(0,8))
temp <- density(myFreq, bw=.05)
temp


tre <- njs(dist(as.matrix(GBS)))
tre$edge.length[tre$edge.length<0]<-0
plot(tre, show.tip=TRUE) # can also use typ="cladogram" or typ="fan"
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
title("Neighbour-joining tree of unfiltered Adneocystis data")
add.scale.bar()

###


my_genind <- vcfR2genind(vcf)
my_genind
head(locNames(my_genind))

my_genind@pop <- factorname

D <- dist(tab(my_genind))
D

tre <- nj(D)
par(xpd=TRUE)
plot(tre, type="unrooted", edge.w=2)
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))

pco1 <- dudi.pco(D) #I chose 3
pco1
dudi.pco(d = D, scannf = FALSE, nf = 7)
s.label(pco1$li*1.0, clab=1, pch=2)
textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li),cex=1.4, new=FALSE, xpd=TRUE)
title("Principal Coordinate Analysis\nunfiltered SNP data")


##


#summary(my_genind@tab)
X <- tab(my_genind, NA.method="mean")

## make PCA
pca1 <- dudi.pca(X,scale=FALSE) #I chose 3
dudi.pca(df = X, scale = FALSE, scannf = FALSE, nf = 3)
temp <- as.integer(pop(my_genind))


# basic plot
plot(pca1$li, cex=3)
library(wordcloud)
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)
## legend the axes by adding loadings
abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1*.5, add.plot=TRUE)
legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7),leg=c("Group A","Group B"), pt.cex=2)





#### 

temp <- inbreeding(my_genind)

class(temp)
head(names(temp))
head(temp[[1]],20)

Fbar <- sapply(temp, mean)

par(mfrow=c(1,1))
hist(Fbar, col="firebrick", main="Average inbreeding in Adenocystis")

which(Fbar>0.5)

F <- inbreeding(my_genind, res.type="function")[which(Fbar>0.5)]
F
#####

ggplot(dpf, aes(x=Sample, y=Depth)) + geom_boxplot(fill=palette) + theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(trans=scales::log2_trans(), expand = c(0,0), breaks=c(1, 10, 100, 1000, 5000),
                     minor_breaks=c(1:10, 2:10*10, 2:10*100, 2:5*1000)) +
  theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6)) +
  theme(panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2))
