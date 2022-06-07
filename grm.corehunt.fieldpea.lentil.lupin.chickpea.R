# select the most diverse lines for Fieldpea, Lentil, lupin, and chickpea
# method:
# load the genotypic data for NVT lines
setwd("/group/grains/pulses/adnan/misc_task/grm_crops")

# Standard method, working on vcf file
library(vcfR)
vcf_len <- read.vcfR("./raw.data/220531_Lentil-CS_SingleHyb.imputationSet.selected.recode.vcf.gz", verbose = F) # variants = 70514; fix_cols = 8; gt_cols = 5083

#########
#***** Object of Class vcfR *****
#  237 samples
# 7 CHROMs
# 8,640 variants
# Object size: 17.5 Mb
# 3.756 percent missing data
#*****        *****         *****
#########

gt.len <- extract.gt(vcf_len, element = "GT", as.numeric = FALSE)
dim(gt.len) # 9640   237
rownames(gt.len) <- paste0("chr", rownames(gt.len))
colnames(gt.len) <- gsub("\\|", "\\.", colnames(gt.len))
colnames(gt.len) <- gsub("-", "\\.", colnames(gt.len))

#gt.len[1:5, 1:5]

# genotype names
library(readxl)
library(janitor)
library(tidyverse)
library(stringr)
library(tidyr)
library(dplyr)

tgt.len <- t(gt.len) #transposing
#dim(tgt.len) # 54 10245
tgtr.len <- cbind(IND = row.names(tgt.len), tgt.len)
#dim(tgtr.len) # 54 10245
row.names(tgtr.len) <- NULL
tgtr.len[1:5, 1:5]

plantIDs.len <- tgtr.len[, 1]
X.len <- tgtr.len[, -1]
X.len[X.len == "0/0"] <- 0
X.len[X.len == "0/1"] <- 1
X.len[X.len == "1/0"] <- 1
X.len[X.len == "1/1"] <- 2
X.len <- apply(X.len, 2, as.numeric)
rownames(X.len) <- plantIDs.len
# print("dim of final X ...")
dim(X.len) # 237 1397
# X.len[1:5,1:5]

# using regrex expression to drop prefix and suffix to identifying the line names  
plantIDs.len <- as_tibble(plantIDs.len) %>%
  separate(value, into = c("chip", "name"), sep = ".CS.") %>%
  separate(name, into = c("chip", "fin"), sep = ".DAV") %>%
  mutate(chip2 = gsub("(.*)\\.\\w+\\d+", "\\1", chip))

len.names <- plantIDs.len$chip2

## calculating minor allele frequency
#retain Xs with MAF of higher than 0.01
maf.len <- (colSums(X.len == 0) + (colSums(X.len == 1) * 0.5)) / nrow(X.len)
# idx <- which(maf>=0.03 & maf<= 0.97)
idx.len <- which(maf.len >= 0.01)
X.len <- X.len[, idx.len]
dim(X.len) # 17 9681

A <- X.len
row.names(A) <- len.names
dim(A)
class(A)
A[1:5, 1:5]
#              chr1:190069 chr1:967538 chr1:992796
# PBA.Jumbo2           0           0           0
# GIA2004L             0           2           2
# GIA2003L             0          NA           2

# calculate distance matrix
a <- dist(A, method = "euclidean")

#            PBA.Jumbo2 GIA2004L GIA2003L
# PBA.Jumbo2    0.00000 52.10892 56.02565
# GIA2004L     52.10892  0.00000 56.64687
# GIA2003L     56.02565 56.64687  0.00000

# why Comparing MDS and PCA
# Mathematically and conceptually, there are close correspondences between MDS and other methods used to
# reduce the dimensionality of complex data, such as Principal components analysis (PCA) and factor analysis.

# PCA is more focused on the dimensions themselves, and seek to maximize explained variance,
# whereas MDS is more focused on relations among the scaled objects.

# MDS projects n-dimensional data points to a (commonly) 2-dimensional space such that
# similar objects in the n-dimensional space will be close together on the two dimensional plot,
# while PCA projects a multidimensional space to the directions of maximum variability
# using covariance/correlation matrix to analyze the correlation between data points and variables.

# perform multi dimensional scaline
fit <- cmdscale(a, eig = TRUE, k = 5) # k is the number of dim
fit # view results
# working with mds results
mds.points <- data.frame(fit[["points"]])
colnames(mds.points) <- c("Dim1", "Dim2", "Dim3", "Dim4", "Dim5")
mds.points$names <- row.names(mds.points)

row.names(mds.points) <- NULL
mds.points$names <- as.factor(mds.points$names)

head(mds.points)

#            names       Dim1       Dim2 Type
# 1      PBA.Jumbo2   4.526275 -14.812524  nvt
# 2        GIA2004L  -3.102435  -1.386384  nvt
# 3        GIA2003L  -6.772216 -20.107828  nvt
# 4   PBA.Kelpie.XT  -6.498693  -0.757372  nvt
# 5 PBA.Hallmark.XT -12.796998   6.482505  nvt
# 6        GIA2001L  -1.296463   2.532126  nvt

options(bitmapType = "cairo")
library(ggplot2)
library(tidyverse)
library(ggrepel)
#dir.create("figures")
# plot lentil lines only
p1 <- mds.points %>%
  ggplot() +
  geom_point( data = mds.points,
              aes(x = Dim1, y = Dim2), size = 2, alpha = 1, col = "blue", alpha = 0.5) +
  geom_text_repel(aes(x = Dim1, y = Dim2, label = ifelse(Dim2 > 10, as.character(names), "")), hjust = 1, vjust = 0) + # present only selective labels
  coord_equal() +
  theme_bw() +
  ggtitle("MDS lentil 237 lines")
p1
ggsave("./figures/lentil_mds.png")

#########################
# Performing Dendrogram analysis
##########################
# perform cluster anlysis
res.hc.complete <- hclust(a)
options(bitmapType = "cairo")

library(dendextend)
dend <- as.dendrogram(res.hc.complete)
# to identify the upper branch of dendrogram
# and further zoom
pdf("./figures/lentil.dend.pdf", width = 40, height = 10)
plot(dend, main = "lentil 237 lines")
dev.off()

#####################
geno.len <- cbind(NAME = row.names(A), A)
dim(geno.len) # 54 10245
row.names(geno.len) <- NULL
geno.len <- cbind(ID = 1:237, geno.len)
geno.len[1:5, 1:5] # final file format used in core hunter

#      ID  NAME                chrchrLentil_Lcu.2RBY.Chr1_424297 chrchrLentil_Lcu.2RBY.Chr1_1996566_O
#[1,] "1" "03.011L.07H4023"   "0"                               "2"                                 
#[2,] "2" "02.427L.07H4028"   "0"                               "0"                                 
#[3,] "3" "03.098L.8.07H4019" "0"                               "2"                                 
#[4,] "4" "03.109L.4.06H4010" "0"                               "2"                                 
#[5,] "5" "03.103L.4.07H4043" "0"                               "2"                                 

write.csv(geno.len, "./clean_files/geno.len.csv", row.names = FALSE)

# install.packages("james.analysis")
# install.packages("rJava")
library(rJava)
library(corehunter)
###
len.data <- coreHunterData(genotypes = genotypes(file = "./clean_files/geno.len.csv", format = "biparental"))

obj <- objective("EN")
len.core <- sampleCore(len.data, obj, size = 55)

len.sel <- len.core[["sel"]]
#saving names and ids of the lines
len.org <- data.frame(cbind(len.data[["ids"]], len.data[["names"]]))
write.csv(len.org, "./clean_files/lentil.name.list.csv", row.names = F)

#saving the names and ids of the core collection
len.sel.name <- len.org[len.org$X1 %in% len.sel, ]
colnames(len.sel.name) <- c("id", "names")
len.sel.name$names <- as.character(len.sel.name$names)
write.csv(len.sel.name, "./clean_files/lentil.core.selection.name.list_EN.csv", row.names = F)

len.sel.mds <- mds.points[mds.points$names %in% len.sel.name$names, ]
head(mds.points)
# plot lentil lines only
p2 <- mds.points %>%
  ggplot() +
  geom_point(data = mds.points, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = len.sel.mds, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  coord_equal() +
  theme_bw() +
  ggtitle("lentil core using Average entry to nearest entry distance (genetic diversity)")
p2
ggsave("./figures/lentil_mds_with_core_EN.png")

################################
obj.AN <- objective("AN")
len.core.AN <- sampleCore(len.data, obj = obj.AN, size = 55)

len.sel.AN <- len.core.AN[["sel"]] # sub-setting the diverse selected lines

len.org <- data.frame(cbind(len.data[["ids"]], len.data[["names"]])) # subsettiname and ids of the orginal lines
# write.csv(len.org, "./clean_files/fieldpea.name.list.csv", row.names = F)

len.sel.name.AN <- len.org[len.org$X1 %in% len.sel.AN, ]
colnames(len.sel.name.AN) <- c("id", "names")
len.sel.name.AN$names <- as.character(len.sel.name.AN$names)
write.csv(len.sel.name.AN, "./clean_files/lentil.core.selection.name.list_AN.csv", row.names = F)

len.sel.mds.AN <- mds.points[mds.points$names %in% len.sel.name.AN$names, ]
# ifelse(Type == "nvt", as.character(names), '')
# plot fieldpea lines only
p2.len.AN <- mds.points %>%
  ggplot() +
  geom_point(data = mds.points, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = len.sel.mds.AN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("fieldpea core using Average accession to nearest entry distance (representativeness)")
p2.len.AN
ggsave("./figures/lentil_mds_with_core_using_AN.png")

##########################################################################
######################### Fieldpea
##########################################################################
# Standard method, working on vcf file
library(vcfR)
vcf_fp <- read.vcfR("./raw.data/220531_Fieldpea-CS_SingleHyb.imputationSet.selected.recode.vcf.gz", verbose = F) # variants = 70514; fix_cols = 8; gt_cols = 5083

gt.fp <- extract.gt(vcf_fp, element = "GT", as.numeric = FALSE)
dim(gt.fp) # 70514   5082
rownames(gt.fp) <- paste0("chr", rownames(gt.fp))
colnames(gt.fp) <- gsub("\\|", "\\.", colnames(gt.fp))
colnames(gt.fp) <- gsub("-", "\\.", colnames(gt.fp))

gt.fp[1:5, 1:5]

tgt.fp <- t(gt.fp)
dim(tgt.fp) # 54 10245
tgtr.fp <- cbind(IND = row.names(tgt.fp), tgt.fp)
dim(tgtr.fp) # 54 10245
row.names(tgtr.fp) <- NULL
tgtr.fp[1:5, 1:5]

plantIDs.fp <- tgtr.fp[, 1]
X.fp <- tgtr.fp[, -1]
X.fp[X.fp == "0/0"] <- 0
X.fp[X.fp == "0/1"] <- 1
X.fp[X.fp == "1/0"] <- 1
X.fp[X.fp == "1/1"] <- 2
X.fp <- apply(X.fp, 2, as.numeric)
rownames(X.fp) <- plantIDs.fp
# print("dim of final X ...")
dim(X.fp)
X.fp[1:5, 1:5]

plantIDs.fp <- as_tibble(plantIDs.fp) %>%
  separate(value, into = c("chip", "name"), sep = ".CS.") %>%
  separate(name, into = c("chip", "fin"), sep = ".DAV") %>%
  mutate(chip2 = gsub("(.*)\\.\\w+\\d+", "\\1", chip))

fp.names <- plantIDs.fp$chip2

## retain Xs with MAF of higher than 0.01
maf.fp <- (colSums(X.fp == 0) + (colSums(X.fp == 1) * 0.5)) / nrow(X.fp)
# idx <- which(maf>=0.03 & maf<= 0.97)
idx.fp <- which(maf.fp >= 0.01)
X.fp <- X.fp[, idx.fp]
dim(X.fp) # 204 937

A.fp <- X.fp
dim(A.fp)
class(A.fp)
A.fp[1:5, 1:5]

row.names(A.fp) <- fp.names

# calculate distance matrix
a.fp <- dist(A.fp, method = "euclidean")
a.fp.1 <- as.matrix(a.fp)
a.fp.1 <- cbind(NAME = row.names(a.fp.1), a.fp.1) # 54 10245
a.fp.1 <- cbind(ID = 1:204, a.fp.1) # 54 10245
row.names(a.fp.1) <- NULL
colnames(a.fp.1)[3:206] <- 1:204
a.fp.1[1:5, 1:5]
write.csv(a.fp.1, "./clean_files/fieldpea_matrix.csv", row.names = F)

# perform multi dimensional scaline
fit.fp <- cmdscale(a.fp, eig = TRUE, k = 5) # k is the number of dim
fit.fp # view results

# working with mds results
mds.points.fp <- data.frame(fit.fp[["points"]])
colnames(mds.points.fp) <- c("Dim1", "Dim2", "Dim3", "Dim4", "Dim5")
mds.points.fp$names <- row.names(mds.points.fp)
row.names(mds.points.fp) <- NULL
mds.points.fp$names <- as.factor(mds.points.fp$names)
head(mds.points.fp)

options(bitmapType = "cairo")
library(ggplot2)
library(tidyverse)
library(ggrepel)
#dir.create("figures")

# plot lentil lines only
p1.fp <- mds.points.fp %>%
  ggplot(aes(x = Dim1, y = Dim2), size = 2, alpha = 1) +
  geom_point(size = 1, col = "blue", alpha = 0.5) +
  geom_text_repel(aes(x = Dim1, y = Dim2, label = ifelse(Dim2 > 25, as.character(names), "")), hjust = 0, vjust = 0) + # present only selective labels
  theme_bw() +
  ggtitle("MDS fieldpea 204 lines")
p1.fp
ggsave("./figures/fieldpea_mds.png")

#########################
# Performing Dendrogram analysis
##########################
# perform cluster anlysis
res.hc.complete.fp <- hclust(a.fp)

options(bitmapType = "cairo")

library(dendextend)
dend.fp <- as.dendrogram(res.hc.complete.fp)
# to identify the upper branch of dendrogram
# and further zoom
pdf("./figures/fieldpea.dend.pdf", width = 40, height = 10)
plot(dend.fp, main = "fieldpea 204 lines")
dev.off()

#####################
geno.fp <- cbind(NAME = row.names(A.fp), A.fp)
dim(geno.fp) # 54 10245
row.names(geno.fp) <- NULL
geno.fp <- cbind(ID = 1:204, geno.fp)
geno.fp[1:5, 1:5]
write.csv(geno.fp, "./clean_files/geno.fp.csv", row.names = FALSE)

# library(rJava)
# library(corehunter)
###
fp.data <- coreHunterData(genotypes = genotypes(file = "./clean_files/geno.fp.csv", format = "biparental"))
#
obj.EN <- objective("EN")
fp.core.EN <- sampleCore(fp.data, obj = obj.EN, size = 55)
fp.sel.EN <- fp.core.EN[["sel"]] # subsetting the diverse selected lines

fp.org <- data.frame(cbind(fp.data[["ids"]], fp.data[["names"]])) # subsettiname and ids of the orginal lines
write.csv(len.org, "./clean_files/fieldpea.name.list.csv", row.names = F)

fp.sel.name.EN <- fp.org[fp.org$X1 %in% fp.sel.EN, ]
colnames(fp.sel.name.EN) <- c("id", "names")
fp.sel.name.EN$names <- as.character(fp.sel.name.EN$names)
write.csv(fp.sel.name.EN, "./clean_files/fieldpea.core.selection.name.list_EN.csv", row.names = F)

fp.sel.mds.EN <- mds.points.fp[mds.points.fp$names %in% fp.sel.name.EN$names, ]
head(mds.points.fp)

# plot fieldpea lines using EN
p2.fp.EN <- mds.points.fp %>%
  ggplot() +
  geom_point(data = mds.points.fp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = fp.sel.mds.EN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("fieldpea core using Average entry to nearest entry distance (genetic diversity)")
p2.fp.EN
ggsave("./figures/fieldpea_mds_with_core_using_EN.png")


#################
# fp.data <- coreHunterData(genotypes = genotypes(file ="./clean_files/geno.fp.csv", format = "biparental"))
obj.AN <- objective("AN")
fp.core.AN <- sampleCore(fp.data, obj = obj.AN, size = 55)

fp.sel.AN <- fp.core.AN[["sel"]] # subsetting the diverse selected lines

fp.org <- data.frame(cbind(fp.data[["ids"]], fp.data[["names"]])) # subsettiname and ids of the orginal lines
# write.csv(len.org, "./clean_files/fieldpea.name.list.csv", row.names = F)

fp.sel.name.AN <- fp.org[fp.org$X1 %in% fp.sel.AN, ]
colnames(fp.sel.name.AN) <- c("id", "names")
fp.sel.name.AN$names <- as.character(fp.sel.name.AN$names)
write.csv(fp.sel.name.AN, "./clean_files/fieldpea.core.selection.name.list_AN.csv", row.names = F)

fp.sel.mds.AN <- mds.points.fp[mds.points.fp$names %in% fp.sel.name.AN$names, ]
head(mds.points.fp)
# plot fieldpea lines only
p2.fp.AN <- mds.points.fp %>%
  ggplot() +
  geom_point(data = mds.points.fp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = fp.sel.mds.AN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("fieldpea core using Average accession to nearest entry distance (representativeness)")
p2.fp.AN
ggsave("./figures/fieldpea_mds_with_core_using_AN.png")

##########################################################################
######################### Chickpea
##########################################################################
# Standard method, working on vcf file
library(vcfR)
vcf_cp <- read.vcfR("./raw.data/220531_Chickpea-CS_SingleHyb.imputationSet.selected.recode.vcf.gz", verbose = F) # variants = 70514; fix_cols = 8; gt_cols = 5083

gt.cp <- extract.gt(vcf_cp, element = "GT", as.numeric = FALSE)
dim(gt.cp) # 70514   5082
rownames(gt.cp) <- paste0("chr", rownames(gt.cp))
colnames(gt.cp) <- gsub("\\|", "\\.", colnames(gt.cp))
colnames(gt.cp) <- gsub("-", "\\.", colnames(gt.cp))

gt.cp[1:5, 1:5]
dim(gt.cp)
colnames(gt.cp)
tgt.cp <- t(gt.cp)
dim(tgt.cp) # 54 10245
tgtr.cp <- cbind(IND = row.names(tgt.cp), tgt.cp)
dim(tgtr.cp) # 54 10245
row.names(tgtr.cp) <- NULL
tgtr.cp[1:5, 1:5]

plantIDs.cp <- tgtr.cp[, 1]
X.cp <- tgtr.cp[, -1]
X.cp[X.cp == "0/0"] <- 0
X.cp[X.cp == "0/1"] <- 1
X.cp[X.cp == "1/0"] <- 1
X.cp[X.cp == "1/1"] <- 2
X.cp <- apply(X.cp, 2, as.numeric)
rownames(X.cp) <- plantIDs.cp
# print("dim of final X ...")
dim(X.cp)
X.cp[1:5, 1:5]

plantIDs.cp <- as_tibble(plantIDs.cp) %>%
  mutate(name = gsub("^.*?CS.", "", value)) %>%
  mutate(name2 = gsub("^.*?DPI.", "", name)) %>%
  mutate(name3 = gsub("(.*)\\.\\w+\\d+\\.\\w+", "\\1", name2)) %>%
  mutate(name4 = gsub(".(A|B|C|D|E|F|G|H)\\d{2}$", "", name3))

cp.names <- plantIDs.cp$name4

# NO MAF applied to  chickpea markers

## retain Xs with MAF of higher than 0.01
# maf.cp <- (colSums(X.cp==0)+(colSums(X.cp==1)*0.5))/nrow(X.cp)
# idx <- which(maf>=0.03 & maf<= 0.97)
# idx.cp <- which(maf.cp>=0.01)
# X.cp <- X.cp[,idx.cp]
# dim(X.cp) # 204 937
X.cp[1:5, 1:5]
A.cp <- X.cp
dim(A.cp)
class(A.cp)
A.cp[1:5, 1:5]

row.names(A.cp) <- cp.names

# calculate distance matrix
a.cp <- dist(A.cp, method = "euclidean")
a.cp.1 <- as.matrix(a.cp)
a.cp.1 <- cbind(NAME = row.names(a.cp.1), a.cp.1) # 54 10245
a.cp.1 <- cbind(ID = 1:391, a.cp.1) # 54 10245
row.names(a.cp.1) <- NULL
colnames(a.cp.1)[3:393] <- 1:391
a.cp.1[1:5, 1:5]
write.csv(a.cp.1, "./clean_files/chickpea_matrix.csv", row.names = F)

# perform multi dimensional scaline
fit.cp <- cmdscale(a.cp, eig = TRUE, k = 5) # k is the number of dim
fit.cp # view results

# working with mds results
mds.points.cp <- data.frame(fit.cp[["points"]])
colnames(mds.points.cp) <- c("Dim1", "Dim2", "Dim3", "Dim4", "Dim5")
mds.points.cp$names <- row.names(mds.points.cp)
row.names(mds.points.cp) <- NULL
mds.points.cp$names <- as.factor(mds.points.cp$names)
head(mds.points.cp)

options(bitmapType = "cairo")
library(ggplot2)
library(tidyverse)
library(ggrepel)
#dir.create("figures")

# plot lentil lines only
p1.cp <- mds.points.cp %>%
  ggplot(aes(x = Dim1, y = Dim2), size = 2, alpha = 1) +
  geom_point(size = 1, col = "blue", alpha = 0.5) +
  geom_text_repel(aes(x = Dim1, y = Dim2, label = ifelse(Dim2 > 25, as.character(names), "")), hjust = 0, vjust = 0) + # present only selective labels
  theme_bw() +
  ggtitle("MDS chickpea 391 lines")
p1.cp
ggsave("./figures/chickpea_mds.png")

#########################
# Performing Dendrogram analysis
##########################
# perform cluster anlysis
res.hc.complete.cp <- hclust(a.cp)

options(bitmapType = "cairo")

library(dendextend)
dend.cp <- as.dendrogram(res.hc.complete.cp)
# to identify the upper branch of dendrogram
# and further zoom
pdf("./figures/chickpea.dend.pdf", width = 40, height = 10)
plot(dend.cp, main = "chickpea 204 lines")
dev.off()

#####################
geno.cp <- cbind(NAME = row.names(A.cp), A.cp)
dim(geno.cp) # 54 10245
row.names(geno.cp) <- NULL
geno.cp <- cbind(ID = 1:391, geno.cp)
geno.cp[1:5, 1:5]
str(geno.cp)
class(geno.cp)
geno.cp[is.na(geno.cp)] <- ""
write.csv(geno.cp, "./clean_files/geno.cp.csv", row.names = FALSE)

# library(rJava)
# library(corehunter)
###
cp.data <- coreHunterData(genotypes = genotypes(file = "./clean_files/geno.cp.csv", format = "biparental"))
#
##
obj.EN <- objective("EN")
cp.core.EN <- sampleCore(cp.data, obj = obj.EN, size = 55)
cp.sel.EN <- cp.core.EN[["sel"]] # subsetting the diverse selected lines

cp.org <- data.frame(cbind(cp.data[["ids"]], cp.data[["names"]])) # subsettiname and ids of the orginal lines
write.csv(cp.org, "./clean_files/chickpea.name.list.csv", row.names = F)

cp.sel.name.EN <- cp.org[cp.org$X1 %in% cp.sel.EN, ]
colnames(cp.sel.name.EN) <- c("id", "names")
cp.sel.name.EN$names <- as.character(cp.sel.name.EN$names)
write.csv(cp.sel.name.EN, "./clean_files/chickpea.core.selection.name.list_EN.csv", row.names = F)

cp.sel.mds.EN <- mds.points.cp[mds.points.cp$names %in% cp.sel.name.EN$names, ]
head(mds.points.cp)

# plot fieldpea lines using EN
p2.cp.EN <- mds.points.cp %>%
  ggplot() +
  geom_point(data = mds.points.cp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = cp.sel.mds.EN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  scale_x_continuous(breaks = c(-200, -150, -100, -50, -25, -10, 0, 10, 20, 30, 40, 50, 150)) +
  scale_y_continuous(limits = c(-300, 210)) +
  theme_bw() +
  ggtitle("chickpea core using Average entry to nearest entry distance (genetic diversity)")
p2.cp.EN
ggsave("./figures/chickpea_mds_with_core_using_EN.png")

#################
# fp.data <- coreHunterData(genotypes = genotypes(file ="./clean_files/geno.fp.csv", format = "biparental"))
obj.AN <- objective("AN")
cp.core.AN <- sampleCore(cp.data, obj = obj.AN, size = 55)

cp.sel.AN <- cp.core.AN[["sel"]] # subsetting the diverse selected lines

cp.org <- data.frame(cbind(cp.data[["ids"]], cp.data[["names"]])) # sub-setting name and ids of the original lines
# write.csv(len.org, "./clean_files/chickpea.name.list.csv", row.names = F)

cp.sel.name.AN <- cp.org[cp.org$X1 %in% cp.sel.AN, ]
colnames(cp.sel.name.AN) <- c("id", "names")
cp.sel.name.AN$names <- as.character(cp.sel.name.AN$names)
write.csv(cp.sel.name.AN, "./clean_files/chickpea.core.selection.name.list_AN.csv", row.names = F)

cp.sel.mds.AN <- mds.points.cp[mds.points.cp$names %in% cp.sel.name.AN$names, ]
head(mds.points.cp)

# plot fieldpea lines only
p2.cp.AN <- mds.points.cp %>%
  ggplot() +
  geom_point(data = mds.points.cp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = cp.sel.mds.AN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("chickpea core using Average accession to nearest entry distance (representativeness)")
p2.cp.AN
ggsave("./figures/chickpea_mds_with_core_using_AN.png")

##########################################################################
######################### LUPIN
##########################################################################
# Standard method, working on vcf file
library(vcfR)
vcf_lp <- read.vcfR("./raw.data/220531_Lupin-CS_SingleHyb.imputationSet.selected.recode.vcf.gz", verbose = F) # variants = 70514; fix_cols = 8; gt_cols = 5083

gt.lp <- extract.gt(vcf_lp, element = "GT", as.numeric = FALSE)
dim(gt.lp) # 70514   5082
rownames(gt.lp) <- paste0("chr", rownames(gt.lp))
colnames(gt.lp) <- gsub("\\|", "\\.", colnames(gt.lp))
colnames(gt.lp) <- gsub("-", "\\.", colnames(gt.lp))

gt.lp[1:5, 1:5]
dim(gt.lp)
colnames(gt.lp)
tgt.lp <- t(gt.lp)
dim(tgt.lp) # 54 10245
tgtr.lp <- cbind(IND = row.names(tgt.lp), tgt.lp)
dim(tgtr.lp) # 54 10245
row.names(tgtr.lp) <- NULL
tgtr.lp[1:5, 1:5]

plantIDs.lp <- tgtr.lp[, 1]
X.lp <- tgtr.lp[, -1]
X.lp[X.lp == "0/0"] <- 0
X.lp[X.lp == "0/1"] <- 1
X.lp[X.lp == "1/0"] <- 1
X.lp[X.lp == "1/1"] <- 2
X.lp <- apply(X.lp, 2, as.numeric)
rownames(X.lp) <- plantIDs.lp
# print("dim of final X ...")
dim(X.lp)
X.lp[1:5, 1:5]

plantIDs.lp <- as_tibble(plantIDs.lp) %>%
  mutate(name = gsub("^.*?CS.", "", value)) %>% # drop the names at .cs.
  mutate(name2 = gsub("(.*)\\.\\w+\\d+\\.\\w+", "\\1", name)) %>%
  mutate(name3 = gsub(".(A|B|C|D|E|F|G|H)\\d{2}$", "", name2))

lp.names <- plantIDs.lp$name3

# NO MAF applied to  chickpea markers

## retain Xs with MAF of higher than 0.01
# maf.lp <- (colSums(X.lp==0)+(colSums(X.lp==1)*0.5))/nrow(X.lp)
# idx <- which(maf>=0.03 & maf<= 0.97)
# idx.lp <- which(maf.lp>=0.01)
# X.lp <- X.lp[,idx.lp]
# dim(X.lp) # 204 937
X.lp[1:5, 1:5]
A.lp <- X.lp
dim(A.lp)
class(A.lp)
A.lp[1:5, 1:5]

row.names(A.lp) <- lp.names

# calculate distance matrix
a.lp <- dist(A.lp, method = "euclidean")
a.lp.1 <- as.matrix(a.lp)
a.lp.1 <- cbind(NAME = row.names(a.lp.1), a.lp.1) # 54 10245
a.lp.1 <- cbind(ID = 1:221, a.lp.1) # 54 10245
row.names(a.lp.1) <- NULL
colnames(a.lp.1)[3:223] <- 1:221
a.lp.1[1:5, 1:5]
write.csv(a.lp.1, "./clean_files/lupin_matrix.csv", row.names = F)

# perform multi dimensional scaline
fit.lp <- cmdscale(a.lp, eig = TRUE, k = 5) # k is the number of dim
fit.lp # view results

# working with mds results
mds.points.lp <- data.frame(fit.lp[["points"]])
colnames(mds.points.lp) <- c("Dim1", "Dim2", "Dim3", "Dim4", "Dim5")
mds.points.lp$names <- row.names(mds.points.lp)
row.names(mds.points.lp) <- NULL
mds.points.lp$names <- as.factor(mds.points.lp$names)
head(mds.points.lp)

options(bitmapType = "cairo")
library(ggplot2)
library(tidyverse)
library(ggrepel)
# dir.create("figures")

# plot lentil lines only
p1.lp <- mds.points.lp %>%
  ggplot(aes(x = Dim1, y = Dim2), size = 2, alpha = 1) +
  geom_point(size = 1, col = "blue", alpha = 0.5) +
  geom_text_repel(aes(x = Dim1, y = Dim2, label = ifelse(Dim2 > 30, as.character(names), "")), hjust = 0, vjust = 0) + # present only selective labels
  theme_bw() +
  ggtitle("MDS lupin 221 lines")
p1.lp
ggsave("./figures/lupin_mds.png")

#########################
# Performing Dendrogram analysis
##########################
# perform cluster anlysis
res.hc.complete.lp <- hclust(a.lp)

options(bitmapType = "cairo")

library(dendextend)
dend.lp <- as.dendrogram(res.hc.complete.lp)
# to identify the upper branch of dendrogram
# and further zoom
pdf("./figures/lupin.dend.pdf", width = 40, height = 10)
plot(dend.lp, main = "lupin 221 lines")
dev.off()

#####################
geno.lp <- cbind(NAME = row.names(A.lp), A.lp)
dim(geno.lp) # 54 10245
row.names(geno.lp) <- NULL
geno.lp <- cbind(ID = 1:221, geno.lp)
#geno.lp[1:5, 1:5]
#str(geno.lp)
#class(geno.lp)
geno.lp[is.na(geno.lp)] <- ""
write.csv(geno.lp, "./clean_files/geno.lp.csv", row.names = FALSE)

# library(rJava)
# library(corehunter)
###
lp.data <- coreHunterData(genotypes = genotypes(file = "./clean_files/geno.lp.csv", format = "biparental"))
##
obj.EN <- objective("EN")
lp.core.EN <- sampleCore(lp.data, obj = obj.EN, size = 55)
lp.sel.EN <- lp.core.EN[["sel"]] # subsetting the diverse selected lines

lp.org <- data.frame(cbind(lp.data[["ids"]], lp.data[["names"]])) # subsettiname and ids of the orginal lines
write.csv(lp.org, "./clean_files/lupin.name.list.csv", row.names = F)

lp.sel.name.EN <- lp.org[lp.org$X1 %in% lp.sel.EN, ]
colnames(lp.sel.name.EN) <- c("id", "names")
lp.sel.name.EN$names <- as.character(lp.sel.name.EN$names)
write.csv(lp.sel.name.EN, "./clean_files/lupin.core.selection.name.list_EN.csv", row.names = F)

lp.sel.mds.EN <- mds.points.lp[mds.points.lp$names %in% lp.sel.name.EN$names, ]
head(mds.points.lp)

# plot lupin lines using EN
p2.lp.EN <- mds.points.lp %>%
  ggplot() +
  geom_point(data = mds.points.lp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = lp.sel.mds.EN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("lupin core using Average entry to nearest entry distance (genetic diversity)")
p2.lp.EN
ggsave("./figures/lupin_mds_with_core_using_EN.png")

#################
# fp.data <- coreHunterData(genotypes = genotypes(file ="./clean_files/geno.fp.csv", format = "biparental"))
obj.AN <- objective("AN")
lp.core.AN <- sampleCore(lp.data, obj = obj.AN, size = 55)

lp.sel.AN <- lp.core.AN[["sel"]] # subsetting the diverse selected lines

lp.org <- data.frame(cbind(lp.data[["ids"]], lp.data[["names"]])) # sub-setting name and ids of the original lines
# write.csv(len.org, "./clean_files/chickpea.name.list.csv", row.names = F)

lp.sel.name.AN <- lp.org[lp.org$X1 %in% lp.sel.AN, ]
colnames(lp.sel.name.AN) <- c("id", "names")
lp.sel.name.AN$names <- as.character(lp.sel.name.AN$names)
write.csv(lp.sel.name.AN, "./clean_files/lupin.core.selection.name.list_AN.csv", row.names = F)

lp.sel.mds.AN <- mds.points.lp[mds.points.lp$names %in% lp.sel.name.AN$names, ]
head(mds.points.lp)

# plot lupin lines only
p2.lp.AN <- mds.points.lp %>%
  ggplot() +
  geom_point(data = mds.points.lp, aes(x = Dim1, y = Dim2), size = 1, col = "blue", alpha = 0.5) +
  geom_point(data = lp.sel.mds.AN, aes(x = Dim1, y = Dim2), shape = 24, fill = "red", size = 2, alpha = 1) +
  theme_bw() +
  ggtitle("lupin core using Average accession to nearest entry distance (representativeness)")
p2.lp.AN
ggsave("./figures/lupin_mds_with_core_using_AN.png")

######################
