library(data.table)
library(tidyverse)
library(ggpubr)

if (!interactive()) {
    # Expected CLI args:
    # 1) .enrich file 2) annotation file
    # 3) output prefix (including path, but excluding file extension)
    args <- commandArgs(trailingOnly = TRUE)
    print(args)
    enrich_file <- args[1]
    annot_file <- args[2]
    output_prefix <- args[3]

    if (length(args) > 3){
        plot_width <- as.numeric(args[4])
        plot_height <- as.numeric(args[5])
    } else{
        plot_width <- 5
        plot_height <- 4
    }

    if (length(args) > 5){
        enrich_cutoff <- as.numeric(args[6])
    } else{
        enrich_cutoff <- 5
    }

    # Load enrichment results
    enrich <- fread(enrich_file)
    colnames(enrich) <- c("Component", "Share", "SE.1", "Expected", "Enrichment", "SE.2", "Z-Stat1", "Z-Stat2")

    # Load annotation file 
    annotations <- fread(annot_file)
    colnames(annotations) <- c("Number","Annotation")

    # Merge them
    enrich <- cbind(enrich,annotations)

    # Remove LDAK_Weightings and Base_Category
    enrich <- enrich %>%
        filter(!grepl("LDAK_Weightings|Base_Category", Annotation))

    # Plot
    enrich %>%
    filter(Enrichment > enrich_cutoff) %>%
    filter(SE.2 != 0) %>%
    mutate(Annotation = factor(Annotation, levels = Annotation[order(Enrichment)])) %>%
    ggplot(aes(x = Annotation, y = Enrichment)) +
    geom_pointrange(aes(ymin = Enrichment - SE.2, ymax = Enrichment + SE.2)) + 
    coord_flip() + 
    theme_light() -> LDAK_enrich

    for (ext in c("pdf", "png")) {
        filename <- paste0(output_prefix, "_LDAK_enrich.", ext)
        if (ext == "png"){
            ggsave(filename, LDAK_enrich, width = plot_width, height = plot_height, units = "in", dpi = 350)
        } else {
            ggsave(filename, LDAK_enrich, width = plot_width, height = plot_height, units = "in")
        }
    }
}
