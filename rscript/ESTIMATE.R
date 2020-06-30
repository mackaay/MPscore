library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)


library(estimate)

OvarianCancerExpr <- system.file("extdata", "sample_input.txt",
                                 package="estimate")
read.table(OvarianCancerExpr)[1:4,1:4]
filterCommonGenes(input.f=OvarianCancerExpr, 
                  output.f="OV_10412genes.gct", 
                  id="GeneSymbol")
common_genes$GeneSymbol

estimateScore(input.ds = "OV_10412genes.gct",
              output.ds="OV_estimate_score.gct", 
              platform="affymetrix")

plotPurity(scores="OV_estimate_score.gct", samples="s516", 
           platform="affymetrix")

scores=read.table("OV_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores


input.f = eSet.file
output.f = "./ESTIMATE/eSet_genes.gct"
id="GeneSymbol"

stopifnot((is.character(input.f) && length(input.f) == 1 && nzchar(input.f)) ||
            (inherits(input.f, "connection") && isOpen(input.f, "r")))
stopifnot((is.character(output.f) && length(output.f) == 1 && nzchar(output.f)))
id <- intersect(common_genes$GeneSymbol,rownames(eSet))   

## Read input data
input.df <- read.table(input.f,
                       header=TRUE,
                       row.names=1,
                       sep="\t", 
                       quote="",
                       stringsAsFactors=FALSE)

merged.df <- merge(common_genes, input.df, by.x="GeneSymbol", by.y="row.names")
rownames(merged.df) <- merged.df$GeneSymbol
merged.df <- merged.df[, -1:-ncol(common_genes)]
print(sprintf("Merged dataset includes %d genes (%d mismatched).",
              nrow(merged.df),
              nrow(common_genes) - nrow(merged.df)))
outputGCT(merged.df, output.f)
