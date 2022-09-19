setwd("~/Documents/Postdoc/Adenocystis_biogeography/DataAnalyses/Adenocystis_Biogeography")

library(vcfR)
library(adegenet)
library(ape)
library(vegan)
library(car)
library(dartR)
library(mmod)
library(reshape2)
library(tidyverse)
library(hierfstat)

### Read in VCF file ###
vcf <- read.vcfR('final.recode.vcf')
head(vcf)
head(getFIX(vcf))


### DEPTH INFO ###

#Extract the depth info
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
dp[1:4,1:6]

dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)

palette=rep_len(c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87"), 99)
# change the end number to however many samples you have

ggplot(dpf, aes(x=Sample, y=Depth)) + geom_boxplot(fill=palette) + theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(trans=scales::log2_trans(), expand = c(0,0), breaks=c(1, 10, 100, 1000, 5000),
                     minor_breaks=c(1:10, 2:10*10, 2:10*100, 2:5*1000)) +
  theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6)) +
  theme(panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2))

### Missingness per sample ###

myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)
myMiss <- data.frame((levels(dpf$Sample)), myMiss)
colnames(myMiss) <- c('Sample', 'Missing')

# change last number to however many samples you have
ggplot(myMiss, aes(x=Sample, y=Missing)) + geom_col(fill=palette) + theme_bw() +
  labs(y='Missingness (%)') +theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=60,hjust=1))+scale_y_continuous(expand = c(0,0))

### Change to a GenLight object ###

GBS<- vcfR2genlight(vcf)
GBS@ind.names

# Best to create a new factorname with whatever you count out
# Always check length!

factorname <- factor(c("BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","Campbell Point","Campbell Point","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Watsons Road North","Watsons Road North","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South"))

factorname <- factor(c("BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek North","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","BullCreek South","Campbell Point","Campbell Point","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kuri Bush South","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Kaka Point","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Lawyer's Head","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks North","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Mitchell Rocks South","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Saint Clair Beach","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Taieri Beach North","Watsons Road North","Watsons Road North","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South","Watsons Road South"))

length(factorname)

GBS@pop <- factorname
GBS@pop
poplevels<-GBS@pop

### Allele frequencies ###

myFreq <- glMean(GBS)
myFreq <- c(myFreq, 1-myFreq)
par(mar = c(5, 5, 5, 5)) # adjust the margins so that they fit on the screen???
hist(myFreq, proba=TRUE, col="#5B88C1", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20,ylim = c(0,8))

temp <- density(myFreq, bw=.05)
temp

par(mfrow=c(1,1))

## Plotting GBS SNPs vs individuals ##
par(mar = c(5, 5, 5, 5))
plot(GBS, col=heat.colors(3), bg="white")

D<-dist(GBS)
class(D)
length(D)
dfD <- as.data.frame(as.matrix(D))
table.paint(dfD, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)
# darker = greater differences
DFd<- t(as.matrix(D))
DFd <- DFd[, ncol(DFd):1]
image(DFd, col = rev(heat.colors(100)))



### Neighbor-joining tree ###

par(mar = c(2, 2, 2, 2))
par(mfrow=c(1,1))

GBS@pop <- factor(GBS@pop, levels = c("Lawyer's Head","Saint Clair Beach","Kuri Bush","Kuri Bush South","Taieri Beach North","Watsons Road North", "Watsons Road South","BullCreek North","BullCreek South","Mitchell Rocks North","Mitchell Rocks South","Kaka Point","Campbell Point"))
levels(GBS@pop)

# basic tree
tree <- njs(dist(as.matrix(GBS)))
par(mar = c(3, 3, 0, 3))
plot(tree, "unrooted", cex=0.45, FALSE, font=4, node.pos=1, edge.width=2, label.offset=1)
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
add.scale.bar()


plot(tree, "phylogram", cex=0.45, FALSE, font=4, node.pos=1, edge.width=2, label.offset=1)
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
axisPhylo()
add.scale.bar()

# with Campbell Point 12 as an outgroup
tree2<-root(tree, out=23)
plot(tree2, "phylogram", cex=0.45, FALSE, font=4, node.pos=1, edge.width=2, label.offset=1)
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
axisPhylo()
add.scale.bar()

# with KP01 as an outgroup
tree2<-root(tree, out=39)
plot(tree2, "phylogram", cex=0.45, FALSE, font=4, node.pos=1, edge.width=2, label.offset=1)
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
axisPhylo()
add.scale.bar()

plot(tree, "fan", cex=0.45, FALSE, font=4, node.pos=1, edge.width=2, label.offset=1)
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
add.scale.bar()

# tree with more pizzaz 
par(mar = c(3, 0, 3, 0))
tre <- njs(dist(as.matrix(GBS)))
tre$edge.length[tre$edge.length<0]<-0
plot(tre, show.tip=TRUE,label.offset=1,font=4,cex=0.45) # can also use typ="cladogram" or typ="fan"
tiplabels(pch=20, cex=1, col=c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")[as.numeric(pop(GBS))])
title("Neighbour-joining tree of Adneocystis data")
axisPhylo()


## PCA w Genlight
GBS
data_pca <- glPca(GBS)

##export list of eigen values and percentage variances for each PC

#eigen values for PCs
eig.val<-data_pca$eig
eig.val
#percentages of variance for each PC
eig.perc <- 100*data_pca$eig/sum(data_pca$eig)
eig.perc
eigen<-data.frame(eig.val,eig.perc)
eigen
#writing file with both
write.csv(eigen,file="eigen-summary.csv",row.names=TRUE,quote=FALSE)

##retrieve the 10 PC scores to be used for plotting later
pca_scores <- as.data.frame(data_pca$scores)
pca_scores$pop <- pop(GBS)
write.table(pca_scores, "10PCs.txt")

##PCA plotting
ggplot(pca_scores, aes(PC1, PC2)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")) +
  theme_classic()

ggplot(pca_scores, aes(PC2, PC3)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")) +
  theme_classic()

ggplot(pca_scores, aes(PC1, PC3)) +
  geom_point(size=3, aes(colour=pop), alpha = 0.7, shape=19) +
  scale_color_manual(values = c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87")) +
  theme_classic()

### GenInd objects ### 

my_genind <- vcfR2genind(vcf)
my_genind@pop <- factorname
my_genind

### Inbreeding ###

inBred <- inbreeding(my_genind)
Fbar <- sapply(inBred, mean)
hist(Fbar, col="firebrick", main="Average inbreeding in Adenocystis")

### Some Stats ###

toto <- summary(my_genind)
names(toto)
par(mar = c(3, 3, 3, 3))
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
par(mar = c(10, 10, 10, 10))
barplot(toto$n.by.pop,ylab="Number of samples",las=3)
par(mar = c(3, 3, 3, 3))

bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")
# observed heterozygosity is lower than mean expected heterozygosity


#summary(my_genind@tab)
X <- tab(my_genind, NA.method="mean")

### F statistics ###
HFS<-pairwise.neifst(genind2hierfstat(my_genind))
FST <- pairwise.fst(HFS, pop=NULL, res.type=c("dist", "matrix"))
fstat(FST)


dapc1 <- dapc(GBS)
scatter(dapc1, scree.da=TRUE,posi.da ="topleft", col = palette)
#Do a DAPC plot; n.pca keeps 100 levels because why not? n.da keeps 13 discriminate things
dapc_pops <- dapc(GBS, n.da = 4, n.pca = 40)
#Graph DAPC plot
scatter(dapc_pops, scree.pca = TRUE, bg="white", pch=20, cstar=0, col=palette, solid=.6,cex=3, clab=0, leg=TRUE, posi.pca="bottomleft")

#Reset the group for realz this time
grp <- find.clusters(GBS, max.n.clust=10, n.pca=40, n.clust =10)

#Check to see what individuals correspond to what group
table(pop(GBS), grp$grp)



#MAKE GENIND OBJECT FROM VCF
x <- vcfR2genind(vcf)
x

#POP NAMES FOR GENIND
pop.names <- factorname
length(pop.names)
pop(x)= pop.names
x

#MAKE HIERFSTAT FROM GENIND
x2 <- genind2hierfstat(x) 

#SUMMARY STATISTICS OVERALL
basicstat <- basic.stats(x2, diploid = TRUE, digits = 2) 
names(basicstat)   

allelic<-allelic.richness(x2,diploid=TRUE)
colnames(allelic$Ar) <- popNames
obs.het<-data.frame(basicstat$Ho)
exp.het<-data.frame(basicstat$Hs)
pop.freq<-data.frame(basicstat$pop.freq)
fis<-data.frame(basicstat$Fis)
perloc<-data.frame(basicstat$perloc)
overall<-data.frame(basicstat$overall)
overall


#SUMMARY STATS PER POP (NA values removed)
obs.het.Mean <- colMeans(obs.het, na.rm = TRUE, dims = 1)
obs.het.Mean
exp.het.Mean <- colMeans(exp.het, na.rm = TRUE, dims = 1)
exp.het.Mean
fis.Mean <- colMeans(fis, na.rm = TRUE, dims = 1)
fis.Mean
allelic.Mean <- colMeans(allelic$Ar, na.rm = TRUE, dims = 1)
allelic.Mean

#CONFIDENCE INTERVALS FOR FIS PER POP
fis.CI <- boot.ppfis(x2)
fis.CI

###### LEA Plots #####

library(LEA)
