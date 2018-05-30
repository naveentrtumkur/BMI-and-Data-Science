## Install packages from Kevin Coombes repository on R-forge (OOMPA)
## see: http://oompa.r-forge.r-project.org/
source("http://silicovore.com/OOMPA/oompaLite.R")
oompaLite()


## install the following package for heatmap.2
install.packages("gplots")

## Need to install packages from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite(c("edgeR"))

## Load the needed packages
library("limma")
library("survival")
library("edgeR")
library("ClassDiscovery")
library("gplots")
library("lattice") 
## library("cvTools")
## library("glmnet")
## library("Hmisc")


## CLINICAL
clin.data <- read.csv("TCGA_HeadNeck_clinical.csv", row.names = 1)
dim(clin.data)

str(clin.data)


## GENE EXPRESSION (mRNA seq)
RNA.counts <- read.csv("TCGA_HeadNeck_mRNAseq_raw_counts.csv", row.names = 1)
dim(RNA.counts)
RNA.counts[1:5, 1:5]


head(rownames(RNA.counts))


## BOXPLOTS OF UNTRANSFORMED DATA
boxplot(RNA.counts[,1:25])

## Not very illuminative since large outliers dominate

## BOXPLOTS OF LOG2 DATA
log2.cnts <- log2(RNA.counts+1)
boxplot(log2.cnts[,1:25])

spca <- SamplePCA(log2.cnts)
plot(spca)

plot(spca, split=clin.data$hpv.status, col = 1:2)
legend("topleft", c("HPV(-)", "HPV(+)"), pch=15, col=1:2)

d.euc <- distanceMatrix(log2.cnts, "euclidean")
## takes a few seconds ...
eucloc <- cmdscale(d.euc)

plot(eucloc[,1], eucloc[,2], pch=15,
     col=as.numeric(clin.data$hpv.status),
     xlab='Coordinate 1', ylab='Coordinate 2',
     main = "MDS plot log2 scale")
legend("topleft", c("HPV(-)", "HPV(+)"), pch=15, col=1:2)


## What about pearson correlation?
d.corr <- distanceMatrix(log2.cnts, "pearson")
corrloc <- cmdscale(d.corr)
plot(corrloc[,1], corrloc[,2], pch=15,
     col=as.numeric(clin.data$hpv.status),
     xlab='Coordinate 1', ylab='Coordinate 2')
legend("topleft", c("HPV(-)", "HPV(+)"), pch=15, col=1:2)


loc3d <- cmdscale(d.corr, k=3)
pc1 <- loc3d[,1]
pc2 <- loc3d[,2]
pc3 <- loc3d[,3]
cloud(pc3 ~ pc1*pc2, pch=15, cex=1.2,
      col=as.numeric(clin.data$hpv.status),
      screen=list(z=55, x=-70),
      perspective=FALSE)

####################################################################
## FILTERING
####################################################################
## NOTE - In LIMMA do filtering PRIOR to scale normalization

## Mean counts per gene
A <- rowMeans(RNA.counts)
median(A)


summary(A)
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##       0.0     108.3    1365.0    6282.0    5005.0 2070000.0
sum(A == 0)
## [1] 61
sum(A < 10)
## [1] 2536
## Filter genes with mean count < 10
idx <- which(A < 10)
RNA.counts <- RNA.counts[-idx, ]
dim(RNA.counts)
## [1] 17966   259
## NOTE - RNA.counts is filtered, but log2.cnts is not yet ... 
log2.cnts <- log2.cnts[-idx, ]
dim(log2.cnts)
## [1] 17966   259


####################################################################
## TMM Normalization
####################################################################

dge <- DGEList(counts=RNA.counts)
dge <- calcNormFactors(dge, method = "TMM")

names(dge)
## [1] "counts"  "samples"
dge$counts[1:5, 1:5]  ## original count data
##       BA.4074 BA.4075 BA.4076 BA.4077 BA.4078
## A1BG     1053     342     411     340     445
## A2BP1      10      27       0      12       2
## A2LD1     552     698     280     475     309
## A2ML1    1492     261   26836   10959    4026
## A2M      8737    4035    4955    8877   17888
head(dge$samples)     ## library sizes and normalization factors

mean.gene <- rowMeans(RNA.counts)
var.gene <- apply(RNA.counts, 1, var)
plot(mean.gene, var.gene)
abline(a=0, b = 1, col=2)

####################################################################
## DESIGN, compare HPV+ vs -, adjusting for batch number
####################################################################

design <- model.matrix(~ 0 + hpv.status + batch.number, data = clin.data)
head(design)
## contrast matrix
contrast.matrix <- makeContrasts(hpv.statusP - hpv.statusN, levels=design)

####################################################################
## VOOM TRANSFORMATION
####################################################################

## Default span value is 0.5, try changing this
v <- voom(dge, design, span = 0.5, plot = TRUE)

## v$E: numeric matrix of normalized expression values on the log2 scale
dim(v$E)
## [1] 17966   259
v$E[1:3,1:5]

## v$weights: numeric matrix of inverse variance weights
dim(v$weights)
## [1] 17966   259
v$weights[1:3,1:5]
##           [,1]      [,2]      [,3]      [,4]      [,5]
## [1,] 0.9976296 0.8681957 0.8684435 0.9227256 0.9171467
## [2,] 0.3167296 0.2782361 0.2783089 0.2722285 0.2927118
## [3,] 1.1566325 1.0006689 1.0009567 1.1036591 1.0583481
## v$design: design matrix
head(v$design)

## v$targets: lib.size
head(v$targets)

## BOXPLOTS BASED ON VOOM TRANSFORMED DATA (v$E)
## v$E = log2-counts per million (logCPM)
boxplot(v$E[,1:25])

####################################################################
## FIT LIMMA MODEL AND EBAYES
####################################################################

fit <- lmFit(v, design)
fitc <- eBayes(contrasts.fit(fit, contrast.matrix))
tt.hpv <- topTable(fitc, coef=1, number = nrow(v), adjust.method = "BY")
## NOTE - BY controls for correlated p-values

## How many with adjusted p-value < 0.05?
sum(tt.hpv$adj.P.Val < 0.05)
## [1] 3218
## 3214, lots to  look at ...

## probably need a stricter p-value threshold
sum(tt.hpv$adj.P.Val < 10^-10)
## [1] 54
## 55

## How many genes pass various p-value filters?  Here use sapply
p.thresholds <- 10^(seq(-2, -10, by = -1))
## p-value thresholds from 10^-2 to 10^-10

sapply(p.thresholds, function(x) sum(tt.hpv$adj.P.Val < x))
## [1] 2116 1204  692  424  272  171  113   84   54
## 2132 1201  699  421  270  170  112   81   55

## Can also incorporate a fold-change threshold (e.g. 2)
sapply(p.thresholds, function(x) sum(tt.hpv$adj.P.Val < x & abs(tt.hpv$logFC) > log(2)))
## [1] 990 667 438 314 224 150 103  78  52
## 998 664 442 309 222 149 102  75  53

####################################################################
## VOLCANO PLOT
####################################################################

## Color-code as red genes with adjusted p-value < 0.01 and
## absolute logFC > log(2)
idx.sig <- which(tt.hpv$adj.P.Val < 0.01 & abs(tt.hpv$logFC) > log(2))
plot(tt.hpv$logFC, -log10(tt.hpv$adj.P.Val), pch = 20,
     xlab = "log fold-change", ylab = "-log10(P-value)",
     main = "Volcano plot of limma results")
points(tt.hpv$logFC[idx.sig], -log10(tt.hpv$adj.P.Val)[idx.sig], pch = 20, col = 2)
