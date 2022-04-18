# David-young

[**YOUNG**]{.ul}

-   The first cell and gene quality control is in [QC_young](src/QC_01_young_jpriller)

-   The normalisation by deconvolution is in [normalise_young](src/normalise_01_young_jpriller)

-   Feature selection and dimensional reduction in [feature_selection_dimred_young](src/feature_selection_dimred_01_young_jpriller)

-   Batch correction in [batchcorrect_young](src/batchcorrect_young_jpriller)

-   Clustering at different resolutions, [clustering_01](src/clustering_01_young_jpriller)

**To Do/ In Progress**

-  First rough annotation in [annotation_01](src/annotation_01_young_jpriller) 

- [ClusterQC](src/cluster_QC_01_young_jpriller) 


-  Stricter cell and gene QC; dimensional reduction and feature selection with the cleaned data[QC_norm_featrueselect_dimred_02](src/QC_norm_featureselect_dimred_02_stop) - **separate cell types**

-  Clustering with different resolutions [clustering_02](src/clustering_02_young_jpriller_stop) - **try to get the BAMs separate from the other immune**


-  Annotation [annotation_02_res0.3](src/annotation_02_young_jpriller_stop)

-  Differential gene expression (DE) between WT and KO and between WT and HETs for each cluster [DE_WT_KO_HET_clusters](src/DE_WT_KO_HET_clusters_young_jpriller_stop) and [DE_WT_HET_microglia](src/DE_WT_HET_microglia_edgeR_stop)

-   Differential abundance of cells between clusters between WT KO and HET:

    -   [DA barplots](src/DA_barplots_young_jpriller_stop)


    -   Using miloR, pairwise comparisons [DA_miloR](src/DA_miloR_young_jpriller_stop)
 

 -   Interactive [Shinyapp](https://annawilliams.shinyapps.io/shinyApp_jpriller) with the analysis until this stage. - **change colours**
 
- remove microglia background in DE?

-  DE using MAST? 
