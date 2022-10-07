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
