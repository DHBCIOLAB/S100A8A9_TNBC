
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")

# Load packages
library(mygene)

## Define pipeline to be run on all datasets
Run <- function(eSet, name)
{
    message('Running analysis on dataset: ', name)

    ## Rank genes by cox proportional hazard regression model "scores"
    genes <- TopCoxGenes(eSet)

    ## Identify the classifier that yields the most optimal
    ## stratification of samples
    classifier <- BestClassifier(eSet, genes)

    ## Draw kaplan meier plots
    PlotKM(Time(eSet), Event(eSet), classifier$groups,
           filename = sprintf('%s_km_plot.pdf', name))

    ## Differential expression analysis
    de <- DiffEq(eSet, classifier$groups, name, poor - good)[[1]]

    message('Done analysing dataset: ', name)

    ## Return results
    list(classifier = classifier, de = de)
}


##
## Main
##
# install_bitbucket("jeevb/survclust")
# install_bitbucket("jeevb/brca")
library(survclust)

## Source datasets
library(brca)

data(ispy1)
data(gse25066)
data(yau)
data(chin)

# Phenotype data is accessed using the 'phenoData' function
phenoData(ispy1)

# Select the samples that are TN with one minor complication. The sample types 'NA' must
# be converted to non-NA values.
tmp.ispy1 <- ispy1$receptor_status == "TN"
tmp.gse25066 <- gse25066$receptor_status == "TN"
tmp.yau <- yau$receptor_status == "TN"
tmp.chin <- chin$receptor_status == "TN"

ind.ispy1 <- union(which(is.na(tmp.ispy1)), which(!tmp.ispy1))
ind.gse25066 <- union(which(is.na(tmp.gse25066)), which(!tmp.gse25066))
ind.yau <- union(which(is.na(tmp.yau)), which(!tmp.yau))
ind.chin <- union(which(is.na(tmp.chin)), which(!tmp.chin))

## Prepare a list of datasets
datasets <- list(
    ispy1     = ExpressionSet2(ispy1[,ind.ispy1],
                               timeCol = 'rfs.t',
                               eventCol = 'rfs.e'),
    gse25066  = ExpressionSet2(gse25066[,ind.gse25066],
                               timeCol = 'drfs_t',
                               eventCol = 'drfs_e'),
    yau       = ExpressionSet2(yau[,ind.yau],
                               timeCol = 't_dmfs',
                               eventCol = 'e_dmfs'),
    chin      = ExpressionSet2(chin[,ind.chin],
                               timeCol = 'disease_time',
                               eventCol = 'disease_binary')
    )

## Run pipeline on datasets
res <- mapply(Run, datasets, names(datasets), SIMPLIFY = FALSE)

# Revised 'res' result file excluding results from ispy1
res.wo.ispy1 <- res[-1]

## Identify genes that are commonly DE between "poor" and "good" prognosis
## groups in at least 3 datasets
common <- CommonDEGs.MC(lapply(res, '[[', 'de'), inAtLeast = 3,
                        iterations = 100000, holm < 0.05)

# Repeat the above analysis after removing ispy1 and groups in at least 2 datasets
common.wo.ispy1 <- CommonDEGs.MC(lapply(res.wo.ispy1, '[[', 'de'), inAtLeast = 2,
                                 iterations = 100000, holm < 0.05)

# Attach gene names to numbers
gene_names <- getGenes(unique(common[[1]]$gene))
scores_with_gene_names <- merge(gene_names, common[[1]], by.x="entrezgene", by.y="gene")
scores_with_gene_names$'_id' <- NULL
scores_with_gene_names$X_score <- NULL
scores_with_gene_names$query <- NULL
scores_with_gene_names$taxid <- NULL

# Repeat attaching gene names to numbers after removing ispy1
gene_names_wo_ispy1 <- getGenes(unique(common.wo.ispy1[[1]]$gene))
scores_with_gene_names_wo_ispy1 <- merge(gene_names_wo_ispy1, common.wo.ispy1[[1]], by.x="entrezgene", by.y="gene")
scores_with_gene_names_wo_ispy1$'_id' <- NULL
scores_with_gene_names_wo_ispy1$X_score <- NULL
scores_with_gene_names_wo_ispy1$query <- NULL
scores_with_gene_names_wo_ispy1$taxid <- NULL

# Order by score and export
scores_with_gene_names <- scores_with_gene_names[with(scores_with_gene_names, order(-score, entrezgene)),]
write.csv(scores_with_gene_names, "output/TN_only_atleast_3_datasets.csv", row.names=FALSE)

# Repeat ordering by scores and exporting the dataset after removing ispy1
scores_with_gene_names_wo_ispy1 <- scores_with_gene_names_wo_ispy1[with(scores_with_gene_names_wo_ispy1, order(-score, entrezgene)),]
write.csv(scores_with_gene_names_wo_ispy1, "output/TN_only_atleast_2_datasets_wo_ispy1.csv", row.names=FALSE)
