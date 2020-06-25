## ----setup, message=FALSE, echo=FALSE------------------------------------
reset_options_scipen <- getOption("scipen")
reset_options_digits <- getOption("digits")
options(scipen = 1, digits = 5)
library(dnapath)
set.seed(12345)

## ------------------------------------------------------------------------
data(meso)
str(meso)

## ----warning=FALSE-------------------------------------------------------
# Run dnapath using the gene expression and group information from meso dataset.
results <- dnapath(meso$gene_expression, 
                   pathway_list = NULL, 
                   groups = meso$groups)
results

## ------------------------------------------------------------------------
plot(results, alpha = 0.05)

## ------------------------------------------------------------------------
data(meso) # Load the gene expression data
data(p53_pathways)

# Run the differntial network analysis.
results <- dnapath(x = meso$gene_expression,
                   pathway_list = p53_pathways,
                   groups = meso$groups)

results

## ------------------------------------------------------------------------
results <- filter_pathways(results, alpha_pathway = 0.2)
results

## ------------------------------------------------------------------------
# The plot layout is stochastic. Setting the RNG seed allows for reproducible plots.
set.seed(123)
plot(results[[1]])

## ------------------------------------------------------------------------
results <- rename_genes(results, to = "symbol", species = "human",
                        dir_save = tempdir())
results[[1]]

## ------------------------------------------------------------------------
set.seed(123) # Reset seed to use same layout as previous plot.
plot(results[[1]])

## ------------------------------------------------------------------------
# Summary table of the edges in pathway 1.
# Note: could instead use the summarize_edges(results[[1]]) function call.
summary(results[[1]], by_gene = FALSE) # Set by_gene to FALSE to obtain edges.

## ------------------------------------------------------------------------
summary(results[[1]], by_gene = FALSE, alpha = 0.05)

## ------------------------------------------------------------------------
plot_pair(results, "FAS", "TP73")

## ----echo=FALSE----------------------------------------------------------
options(scipen = reset_options_scipen, digits = reset_options_digits)

