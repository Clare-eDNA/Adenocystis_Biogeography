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

vcf <- read.vcfR("populations.snps.FirstAttempt.vcf", verbose = FALSE)
plot(vcf)
