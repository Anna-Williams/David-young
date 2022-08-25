
genotypehetwt <- sce$genotype %in% c("HET","WT")
milo_het <- Milo(sce[,genotypehetwt])
colData(milo_het)[["genotype"]] <- droplevels(colData(milo_het)[["genotype"]])

### Construct KNN graph
#
#We need to add the KNN graph to the Milo object. This is stored in the `graph` slot, in [`igraph`](https://igraph.org/r/) format. The `miloR` package includes functionality to build and store the graph from the PCA dimensions stored in the `reducedDim` slot. In this case, we specify that we want to build the graph from the MNN corrected PCA dimensions.
#
#For graph building you need to define a few parameters:
#  
#  - `d`: the number of reduced dimensions to use for KNN refinement. We recommend using the same $d$ used for KNN graph building.In our corrected dataset we have 23 dimensions.
#- `k`: this  affects the power of DA testing, since we need to have enough cells from each sample represented in a neighbourhood to estimate the variance between replicates. On the other side, increasing $k$ too much might lead to over-smoothing. We suggest to start by using the same value for $k$ used for KNN graph building for clustering and UMAP visualization. We will later use some heuristics to evaluate whether the value of $k$ should be increased.
#
#```{r}
milo_het <- buildGraph(milo_het, k = 30, d = 23, reduced.dim = "corrected")
#```
#
### Defining representative neighbourhoods on the KNN graph
#
#We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don't test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by [Gut et al. 2015](https://www.nature.com/articles/nmeth.3545). 
#
#As well as $d$ and $k$, for sampling we need to define a few additional parameters:
#
#- `prop`: the proportion of cells to randomly sample to start with. We suggest using `prop=0.1` for datasets of less than 30k cells. 
#<!-- For bigger datasets using `prop=0.05` should be sufficient (and makes computation faster). -->
#- `refined`: indicates whether you want to use the sampling refinement algorithm, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using `random` instead, is if you have batch corrected your data with a graph based correction algorithm, such as [BBKNN](https://github.com/Teichlab/bbknn), but the results of DA testing will be suboptimal.
#
#```{r}
milo_het <- makeNhoods(milo_het, prop = 0.1, k = 30, d=23, refined = TRUE, reduced_dims = "corrected")
#```
#
#Once we have defined neighbourhoods, we plot the distribution of neighbourhood sizes (i.e. how many cells form each neighbourhood) to evaluate whether the value of $k$ used for graph building was appropriate. We can check this out using the `plotNhoodSizeHist` function. 
#
#As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples. (5x12=60)
#

plotNhoodSizeHist(milo_het)
 
#
#
### Counting cells in neighbourhoods
#
#_milo_het_ leverages the variation in cell numbers between replicates for the same experimental condition to test for differential abundance. Therefore we have to count how many cells from each sample are in each neighbourhood. We need to use the cell metadata and specify which column contains the sample information.
#
#```{r milo_het}
milo_het <- countCells(milo_het, meta.data = as.data.frame(colData(milo_het)), sample="Sample")
#```{r}
head(nhoodCounts(milo_het))
#```{r}
design_het <- data.frame(colData(milo_het))[,c("Sample", "batch", "genotype")]
design_het$batch <- as.factor(design_het$batch)
design_het <- distinct(design_het)
rownames(design_het) <- design_het$Sample
#contrastKO <- c("genotypeKO - genotypeWT") # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>
# we need to use the ~ 0 + Variable expression here so that we have all of the levels of our variable as separate columns in our model matrix
contrastHET <- c("genotypeHET - genotypeWT")
#```{r long-step}
milo_het <- calcNhoodDistance(milo_het, d=23, reduced.dim = "corrected")
#```
da_results_het <- testNhoods(milo_het, design = ~ 0 + genotype + batch, design.df = design_het, model.contrasts = "genotypeHET - genotypeWT", reduced.dim = "corrected")
