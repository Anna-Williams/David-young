#SCRIPT FOR PLOT WITH FC ON Y AXIS, GROUPS IN X AXIS

## set up ---

library(here) #reproducible paths
library(tidyverse) # manipulate dfs and ggplots

# functions
plot_FC_bycluster_ambient <- function(de_results,
                                    ambient = TRUE,
                                    ambient_threshold = 0.1,
                                    FDR_threshold = 0.05,
                                    cols = c("darkred", "darkgreen", "#999999" ),
                                    data = FALSE){
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(de_results, 
                           function(x) {
                             na.omit(as.data.frame(x)) %>% 
                               rownames_to_column("gene")
                           }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>% 
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold, 
                                yes = ifelse(logFC > 0, "upregulated", "downregulated"),
                                no = "not-significant")) 
  
  if(ambient==TRUE){
    # filter the genes that are more than 10% ambient (or other threshold)
    de_results_df <- de_results_df %>% 
      filter(minAmbient < ambient_threshold)
  }
  
  if(data == FALSE){
    ## plot ---
    de_results_df %>% 
      # order the levels to display in correct order in plot
      mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>% 
      # mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>% 
      # sort so the significant values are plotted on top of the non significants
      arrange(desc(Significant)) %>% 
    ggplot(mapping = aes(x=cluster, y=logFC, color=Significant)) + 
      geom_jitter() + 
      scale_colour_manual(values = cols) +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.title = element_blank()) +
      ylab("logFC fire-mice ko to wt") +
      xlab(element_blank())
  }else{
    de_results_df
  }
}

get_ambient_genes <- function(de_results,
                              ambient_threshold = 0.1,
                              FDR_threshold = 0.05){
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(de_results, 
                           function(x) {
                             na.omit(as.data.frame(x)) %>% 
                               rownames_to_column("gene")
                           }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>% 
    filter(minAmbient > ambient_threshold) %>% 
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold, 
                                yes = ifelse(logFC > 0, "upregulated", "downregulated"),
                                no = "not-significant")) %>% 
    # mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>% 
    # sort so the significant values are plotted on top of the non significants
    arrange(desc(Significant))
}
# project
project <- "fire-mice"


## import and transform data ---

# this is the output from edgeR for sc
# a list with the DE results, each element named as one of the clusters and 
# contain a DFrame with the DE for that cluster
de_results <- readRDS(here("processed", project, "DE_results_ambient_edgeR_KO.RDS"))



plot_FC_bycluster_ambient(de_results)
ggsave(here("outs", project, "DE_edgeR", "plots", "FC_clusters_minAmbient.pdf"), 
       height = 7, width = 10)

de_results_df <- plot_FC_cluster_ambient(de_results, data = TRUE)
write.csv(filter(de_results_df, Significant != "not-significant"),
          here("outs", project, "DE_edgeR", "de_results_KO", "de_results_KO_combined_significant_ambientremoved.csv" ),
          row.names = FALSE)


## save csv with deleted genes
# merge all the dfs from the list in a single big dataframe
de_results_df_ambient <- get_ambient_genes(de_results)

write.csv(filter(de_results_df_ambient, Significant != "not-significant"),
          here("outs", project, "DE_edgeR", "de_results_KO", "de_results_KO_combined_significant_onlyambient.csv" ),
          row.names = FALSE)
  