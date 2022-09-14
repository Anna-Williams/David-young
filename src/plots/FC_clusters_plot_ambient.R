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
                                      cols = c("darkred", "darkgreen", "#999999"),
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
      no = "not-significant"
    ))

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
      ggplot(mapping = aes(x = cluster, y = logFC, color = Significant)) +
      geom_jitter() +
      scale_colour_manual(values = cols) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank()
      ) +
      ylab("logFC fire-mice ko to wt") +
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
    # mutate(cluster = forcats::fct_relevel(cluster, c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "OligoAstro", "Oligo_1", "Oligo_2", "OPCs", "mNeurons", "iNeurons & NRPs", "Lymphocytes", "BAMs", "DCs",  "Endothelial", "fEndothelia", "Mural_cells", "ChP_epithelial"))) %>%
    # sort so the significant values are plotted on top of the non significants
    arrange(desc(Significant))
}


# project
project <- "fire-mice"


# loop through the 3 genotype comparisons

for (gnt in c("KO", "HET", "KOvsHET")) {


  # load output from edgeR for each genoype

  de_results_gnt <- readRDS(here("processed", project, paste0("DE_results_ambient_edgeR_", gnt, ".RDS")))

  # the output from edgeR for sc is
  # a list with the DE results, each element named as one of the clusters and
  # contain a DFrame with the DE for that cluster

  plot_FC_bycluster(de_results_gnt, ambient=TRUE, ambient_threshold = 0.25)
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