# SCRIPT FOR PLOT WITH FC ON Y AXIS, GROUPS IN X AXIS

## set up ---

library(here) # reproducible paths
library(tidyverse) # manipulate dfs and ggplots

# functions
plot_FC_bycluster <- function(de_results,
                              ambient = FALSE,
                              ambient_threshold = 0.1,
                              FDR_threshold = 0.05,
                              logFC_threshold = 0,
                              data = FALSE) {
  
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
    function(x) {
      na.omit(as.data.frame(x)) %>%
        rownames_to_column("gene")
    }
  )
  # merge all the dfs from the list in a single big dataframe
  de_results_df <- bind_rows(de_results_dfs, .id = "cluster") %>%
    # add column that indicates if the result is significant
    mutate(Significant = ifelse(FDR < FDR_threshold,
                                yes = ifelse(logFC > logFC_threshold, "upregulated", "downregulated"),
                                no = "not-significant")) %>% 
    mutate(logFC_colour = ifelse(Significant == "not-significant", NA, logFC))  %>%
    mutate(cluster = forcats::fct_relevel(cluster, unique(cluster)))
    
  
  if (ambient == TRUE) {
    # filter the genes that are more than 10% ambient (or other threshold)
    de_results_df <- de_results_df %>%
      filter(minAmbient < ambient_threshold)
  }
  
  if (data == FALSE) {
    ## plot ---
    de_results_df %>%
      # order the levels to display in correct order in plot
      mutate(Significant = forcats::fct_relevel(Significant, c("upregulated", "downregulated", "not-significant"))) %>%
      # mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>%
      # sort so the significant values are plotted on top of the non significants
      arrange(desc(Significant)) %>%
      ggplot(mapping = aes(x = cluster, y = logFC, color = logFC_colour)) +
      geom_jitter() +
      scale_colour_gradient2(high = scales::muted("red"), mid = "white", low = scales::muted("blue")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()
      ) +
      ylab("Log Fold Change") +
      xlab(element_blank())
  } else {
    de_results_df
  }
}

get_ambient_genes <- function(de_results,
                              ambient_threshold = 0.1,
                              FDR_threshold = 0.05) {
  # transform all the DFrame to a standart df format
  de_results_dfs <- lapply(
    de_results,
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
                                no = "not-significant"
    )) %>%
    mutate(cluster = forcats::fct_relevel(cluster, names(de_results))) %>% 
    # sort so the significant values are plotted on top of the non significants
    arrange(desc(Significant))
}


# project
project <- "fire-mice"


# loop through the 2 genotype comparisons

for (gnt in c("KO", "KOvsHET")) {
  
  
  # load output from edgeR for each genoype
  
  de_results_gnt <- readRDS(here("processed", project, paste0("DE_results_ambient_edgeR_", gnt, ".RDS")))
  
  # the output from edgeR for sc is
  # a list with the DE results, each element named as one of the clusters and
  # contain a DFrame with the DE for that cluster
  
  # sort the order from the list
  de_results_gnt <-
  de_results_gnt[c("Astro_1", "Astro_2", "Astro_Oligo", "OPC", "pOPC", "iOligo", "mOligo_1", "mOligo_2", "mOligo_3", "Endothelial", "Mural_cells", "iNeurons_&_NRPs", "mNeuron_ex", "mNeuron_in")]
  
  plot_FC_bycluster(de_results_gnt, ambient=TRUE, ambient_threshold = 0.25) + theme(legend.direction="horizontal")
  ggsave(here("outs", project, "DE_edgeR", "plots", paste0("FC_clusters_", gnt, ".pdf")),
         height = 7, width = 10
  )
  
  de_results_df <- plot_FC_bycluster(de_results_gnt, ambient=TRUE, data = TRUE)
  write.csv(filter(de_results_df, Significant != "not-significant"),
            here("outs", project, "DE_edgeR", paste0("de_results_", gnt), paste0("de_results_", gnt, "_combined_significant_ambientremoved.csv")),
            row.names = FALSE
  )
  
  ## save csv with deleted genes
  de_results_df_ambient <- get_ambient_genes(de_results_gnt)
  
  write.csv(filter(de_results_df_ambient, Significant != "not-significant"),
            here("outs", project, "DE_edgeR", paste0("de_results_", gnt), paste0("de_results_", gnt, "_combined_significant_onlyambient.csv")),
            row.names = FALSE
  )
}

## for the HETs I'll add the comparisons made aside. BAMs and Microglia. 
de_results_het <- readRDS(here("processed", project, paste0("DE_results_ambient_edgeR_", "HET", ".RDS")))

de_results_het[["Microglia"]] <- read.csv(here("outs", project, "DE_edgeR","de_results_HET", "de_resutls_HET_microglia.csv"), 
                                          row.names = 1) %>%
                                  mutate(minAmbient = 0)
de_results_het[["BAMs"]] <- read.csv(here("outs", project, "DE_edgeR","de_results_HET", "de_resutls_HET_BAMs.csv"), 
                                     row.names = 1) %>%
                            mutate(minAmbient = 0)
# sort the order from the list, for plot
de_results_het <-
  de_results_het[c("Microglia", "BAMs", "Astro_1", "Astro_2", "Astro_Oligo", "OPC", "pOPC", "iOligo", "mOligo_1", "mOligo_2", "mOligo_3", "Endothelial", "Mural_cells", "iNeurons_&_NRPs", "mNeuron_ex", "mNeuron_in")]

plot_FC_bycluster(de_results_het, ambient=TRUE, ambient_threshold = 0.25)
ggsave(here("outs", project, "DE_edgeR", "plots", paste0("FC_clusters_", "HET", ".pdf")),
       height = 7, width = 11
)

